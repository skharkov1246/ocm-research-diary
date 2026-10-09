# Раунд 2. Вариант «газ (не энергия) генераторам покупателя»: подготовка газа на факеле, компрессия,
# труба 8-33 км до завода, PRMS у покупателя; покупатель сам ставит газовые двигатели или dual-fuel комплекты.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном (opt/base/pess).
# Поправки аудита NaCN: капитал x2.4/2.6/2.8 к пакетам оборудования; фикс. опекс >= скрин x1.7/2.0/2.4;
# персонал с экспатом; ТОиР 3.5-5% на оборудование; выход по нижней границе; спрос /1.3-2.7 до контракта.
# Второй год: стройка 1-2 г, разгон, загрузка факела 0.78 (VIIRS), падение дебита 8%/г, оборотный капитал,
# девальвация найры 8-16%/г на лаге пересчёта и дебиторке (цена в $, оплата в найре по NFEM), налог 34%.
# Мощность: 4.40 МВт брутто на 1 MMscf/д (исправленная, LHV-КПД 38%) [О, round2/deal/deal_r2.py].
# Запуск: python3 g2g.py -> g2g_out.json и g2g_stdout.txt рядом.
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
TAX = 0.34            # CIT 30% + development levy 4% (NTA 2025) [В]
FX = 1330.0           # ₦/$ 10.2026 [В]
YEARS = 10            # лет эксплуатации = срок GSA с покупателем [Д]
HHV = 1.05            # MMBtu/Mscf [Д 1.03-1.10]
MW_PER = 4.40         # МВт брутто на 1 MMscf/д [О]

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.99, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

# ---------------- труба: стоимость за км ----------------
# Parker N. (2004) UCD-ITS-RR-04-35, регрессия по 893 проектам OGJ 1991-2003, доллары 2000 г. [П]:
#   Construction Cost = [674 D^2 + 11 754 D + 234 085] * L(миль) + 405 000   (D в дюймах; материалы+работы+прочее+ROW)
# Пересчёт в 2026 г.: x1.88 (CPI США 2000->2026, нижняя граница) ... x3.0 (опережающий рост стоимости
# трубопроводного строительства по OGJ) [Д]; нигерийский коэффициент 0.8-1.3 (дешевле труд, дороже импортная труба,
# переходы рек/болот, община, охрана) [Д]. Нигерийской цены малого диаметра за км НЕ НАЙДЕНО (поиск 09.10.2026):
# есть только CSJ 2015 «$0.8-2 млн/км для больших диаметров» [В] и AKK 40" $4.3-5.4 млн/км [В].
def pipe_cost(D, km, esc, ng):
    mi = km/1.609344
    return ((674*D**2 + 11754*D + 234085)*mi + 405000)*esc*ng/1e6

def weymouth_D(q_mmscfd, km, p1=35.0, p2=10.0):
    """минимальный стандартный диаметр (дюймы, внутр.) под пиковый расход; Weymouth, E=0.92, G=0.65, 30°C [О]"""
    mi = km/1.609344; P1, P2 = p1*14.504, p2*14.504
    for D_nom, D in ((3, 3.068), (4, 4.026), (6, 6.065), (8, 7.981), (10, 10.02)):
        q = 433.5*0.92*(520/14.73)*((P1**2-P2**2)/(0.65*546*mi*0.9))**0.5*D**2.667/1e6
        if q >= q_mmscfd: return D_nom, round(q, 2)
    return 12, None

def hp_per_mmscfd(p_suc, p_dis, stages):
    R = (p_dis/p_suc)**(1/stages); F = {2: 1.08, 3: 1.10}[stages]
    return 22*R*stages*F            # GPSA Engineering Data Book, гл.13, BHP/MMscfd [В]

SC = {
 'opt':  dict(lf=0.85, decl=0.05, up=0.97, shrink=0.045, price=8.0, pdrift=0.025, dev_fx=0.08, lag_q=0.25, dso=30, coll=0.98,
              gas_in=0.25+0.10, top=0.0, usd_hp=1000, cm=2.4, dehy=0.30, hcdp=0.20, prms=0.20, esc=1.88, ng=0.8,
              route=1.2, dev=0.6, build=1, ramp=0.8, staff=0.45, sec=0.25, surv=0.006, comm=0.10, ovh=0.15,
              tor_eq=0.035, tor_pipe=0.015, ins=0.010, scr=1.7, demand=1/1.3, p_suc=3.0),
 'base': dict(lf=0.78, decl=0.08, up=0.92, shrink=0.055, price=6.0, pdrift=0.02, dev_fx=0.12, lag_q=0.25, dso=45, coll=0.95,
              gas_in=0.60*HHV+0.15, top=0.7, usd_hp=1250, cm=2.6, dehy=0.35, hcdp=0.30, prms=0.30, esc=2.35, ng=1.0,
              route=1.3, dev=0.8, build=1, ramp=0.7, staff=0.55, sec=0.35, surv=0.010, comm=0.15, ovh=0.20,
              tor_eq=0.0425, tor_pipe=0.020, ins=0.015, scr=2.0, demand=1/2.0, p_suc=2.0),
 'pess': dict(lf=0.60, decl=0.12, up=0.85, shrink=0.07, price=4.5, pdrift=0.0, dev_fx=0.16, lag_q=0.25, dso=90, coll=0.85,
              gas_in=1.00*HHV+0.30, top=0.9, usd_hp=1500, cm=2.8, dehy=0.45, hcdp=0.40, prms=0.45, esc=3.0, ng=1.3,
              route=1.5, dev=1.0, build=2, ramp=0.5, staff=0.65, sec=0.50, surv=0.015, comm=0.30, ovh=0.30,
              tor_eq=0.05, tor_pipe=0.025, ins=0.020, scr=2.4, demand=1/2.7, p_suc=1.5),
}
# price: цена газа у ворот покупателя, $/MMBtu. Потолок -- доставленный CNG/СПГ на адрес покупателя:
#   Greenville СПГ юг ₦230/scm (04.2025) ~ $4.7-4.9; CNG ₦380-520/scm ~ $8.0-10.75 [В, фазы 1-2] -> база $6.0 [Д 4.5-8.0].
# gas_in: газ с факела $/Mscf «всё включено» (цена NGFCP/GCA + плата оператору) [Д, как в модели сделки].
# top: take-or-pay по газу на долю DCQ (только законтрактованный объём, урок раунда 1).
# NMDPRA [П, S.I. No.40 of 2024]: Wholesale Gas Supply licence $0.03/Mscf; MDGIF 0.5% оптовой цены (PIA s.52) [П].
LIC_MSCF = 0.03; MDGIF = 0.005

def nmdpra_fees(Qdes, km, cluster):
    """разовые и годовые сборы NMDPRA, $ млн [П, Fees Regulations 2024: 2-е, 7-е, 9-е приложения]"""
    big = Qdes > 6
    proc = (50e3+2e3 + 20e3+2e3 + 100e3) if big else (10e3+2e3 + 10e3+2e3 + 10e3)      # LTE+LTC+LTO переработки
    pipe = max(500, 50*km) + 1000 + 500 + 1500 + 500*km                                  # PTS, LTE, LTC, PLL, ATC-сбор 500/км
    dist = (10e3 + 100e3) if cluster else 0.0                                            # EOI + лицензия на распределение
    once = (proc + pipe + dist + 100e3/FX*1)/1e6                                         # + Permit to Access Flare Gas ₦100k
    annual = (min(10e3, 500*km) + 1.5e3*math.ceil(Qdes/6) + (50e3 if cluster else 0))/1e6
    return once, annual

def design(Q, dem_mmscfd, n_buyers, km_route, p):
    """Qdes: проектный вход на компрессию = min(доступно в год 1, спрос/(1-усушка)) x 1.3 (пик) [Д]"""
    avail1 = Q*p['lf']*p['up']
    q_in = min(avail1, dem_mmscfd/(1-p['shrink']))
    Qdes = max(0.3, q_in*1.3)
    D, qcap = weymouth_D(Qdes, km_route)
    hpq = hp_per_mmscfd(p['p_suc'], 35.0, 3 if 35.0/p['p_suc'] > 12 else 2)
    hp = hpq*Qdes*1.5                                                   # 3x50% (N+1) [Д]
    pkg = dict(inlet_meter=0.25+0.03*Qdes, compr=hp*p['usd_hp']/1e6, teg=p['dehy']*Qdes**0.6,
               hcdp=p['hcdp']*Qdes**0.6, prms=p['prms']*n_buyers, scada=0.15)
    eq = sum(pkg.values())*p['cm']                                      # аудит x2.4-2.8 к пакету
    pipe = pipe_cost(D, km_route, p['esc'], p['ng'])
    cluster = n_buyers > 1
    fee_once, fee_yr = nmdpra_fees(Qdes, km_route, cluster)
    inst = eq + pipe + p['dev'] + fee_once
    idc = inst*0.14*p['build']/2
    screen = sum(pkg.values()) + km_route*0.15 + 0.3                     # промоутерский скрин: пакет x1, труба $150k/км [Д]
    return dict(Qdes=round(Qdes, 2), D=D, qcap=qcap, hp=round(hp), pkg=round(sum(pkg.values()), 2), eq=round(eq, 2),
                pipe=round(pipe, 2), pipe_per_km=round(pipe/km_route, 3), fees=round(fee_once, 3), fee_yr=fee_yr,
                total=round(inst+idc, 2), screen=round(screen, 2))

def run(Q, dem, n_buyers, km_route, sc='base', ov=None, contract=True, edti=False, naira_price=False,
        conv_lag=0, years=YEARS):
    """Q: факел, MMscf/д (VIIRS 2025); dem: спрос покупателя(ей) на товарный газ, MMscf/д в среднем;
    conv_lag: лет задержки выхода покупателя на полный отбор (поставка двигателей, конверсия) [Д]."""
    p = dict(SC[sc]); p.update(ov or {})
    d = design(Q, dem, n_buyers, km_route, p); C = d['total']
    demand = dem*(1.0 if contract else p['demand'])
    eq_inst = d['eq']; pipe = d['pipe']
    sc_mult = 1.0 if d['Qdes'] <= 2 else (1.2 if d['Qdes'] <= 7 else 1.5)
    fixed_bu = ((p['staff']+p['ovh'])*sc_mult + p['sec'] + p['surv']*km_route + p['comm'] + d['fee_yr']
                + eq_inst*p['tor_eq'] + pipe*p['tor_pipe'] + C*p['ins'])
    fixed_scr = 0.20*sc_mult + 0.10 + d['screen']*(0.02+0.005)
    fixed = max(fixed_bu, fixed_scr*p['scr'])
    b = p['build']; cf = [-C/b]*b
    fxh = 1 - p['dev_fx']*(p['lag_q']/2 + p['dso']/365)               # потеря на лаге пересчёта ₦/$ и дебиторке
    prev = 0.0; rows = []
    for t in range(1, years+1):
        T = t + b                     # годы от 2025 (VIIRS): FID ~2026-27, первый газ через b лет [Д]
        avail = Q*p['lf']*(1-p['decl'])**(T)*p['up']                   # MMscf/д на входе (снижение с 2025 г.)
        ramp = p['ramp'] if t == 1 else 1.0
        if conv_lag and t <= conv_lag: ramp = min(ramp, 0.5)
        sales = min(avail*(1-p['shrink']), demand)*ramp                  # товарный газ, MMscf/д
        take = sales/(1-p['shrink'])
        paid = max(take, p['top']*min(avail, demand/(1-p['shrink'])))     # ToP только на DCQ, не выше доступного
        price = p['price']*((1/(1+p['dev_fx'])**t) if naira_price else (1+p['pdrift'])**(t-1)*fxh)   # pdrift: эскалатор $-цены (CPI США) [Д 0-2.5%/г]
        rev_inv = sales*1000*365*HHV*price/1e6
        rev = rev_inv*p['coll']
        gas = paid*1000*365*(p['gas_in']*1.02**(T-1))/1e6 + sales*1000*365*LIC_MSCF/1e6 + MDGIF*rev_inv
        fx_t = fixed*1.03**(T-1)
        e = rev - gas - fx_t
        dep = C/10
        tx = max(0.0, (e-dep)*TAX)
        ed = (min(tx, 0.05*C) if t <= 5 else 0.0) if edti else 0.0
        wc = p['dso']/365*(rev-prev); prev = rev
        c = e - tx + ed - wc + (p['dso']/365*rev if t == years else 0)
        cf.append(c)
        rows.append(dict(t=t, avail=round(avail, 2), sales=round(sales, 2), rev=round(rev, 2), gas=round(gas, 2),
                         fixed=round(fx_t, 2), ebitda=round(e, 2)))
    i = irr(cf)
    acc = 0; pb = None
    for k, c in enumerate(cf):
        if pb is None and k > 0 and acc < 0 and acc+c >= 0: pb = round(k-1+(-acc)/c - (b-1), 1)
        acc += c
    return dict(Q=Q, dem=dem, n=n_buyers, km=km_route, sc=sc, contract=contract, design=d, fixed=round(fixed, 2),
                opex_vs_screen=round(fixed_bu/fixed_scr, 2), y1=rows[0], y2=rows[1], y5=rows[4],
                npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2), irr=None if i is None else round(i*100, 1),
                payback=pb)

def solve(fn, lo, hi, it=70):
    flo, fhi = fn(lo), fn(hi)
    if (flo > 0) == (fhi > 0): return None
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 2)

# ---------------- спрос покупателя ----------------
E_ANCHOR = 61.9e3       # МВт·ч/г: якорь 11 МВт нетто x 0.72 x 0.92 (как в «энергии якорю») [О фаз 1-2]
def dem_engines(mwh):   # газовые двигатели покупателя, КПД как у нашей станции (4.40 МВт/MMscf/д, собств. нужды 4%)
    return mwh/(MW_PER*0.96*8760)
def dem_dual(mwh, subst=0.55, eff_d=0.36, pen=1.03):
    # dual-fuel: замещение 45-60% в поле (Cat DGB 60%, ComAp 45%), до 70% по каталогу [В]; КПД дизеля 36% LHV [Д]
    return mwh/(MW_PER*0.96*8760)*(0.38/eff_d)*subst*pen
D_ENG = dem_engines(E_ANCHOR); D_DUAL = dem_dual(E_ANCHOR)

# ---------------- экономика покупателя ----------------
def buyer(price, capex_kw=1583, om=0.018, crf_r=0.15, n=10, lf=0.72*0.92):
    """LCOE своей газовой генерации покупателя, $/кВт·ч; капитал $1 500-1 667/кВт [В, фазы 1-2], O&M $0.015-0.02 [Д]"""
    crf = crf_r/(1-(1+crf_r)**-n)
    cap = capex_kw*crf/(8760*lf)
    hr = 3.6/0.38/1.055/0.903/1000*1.0      # MMBtu HHV на кВт·ч: 3.6 МДж/0.38 -> LHV, /0.903 -> HHV [О]
    return cap + om + price*hr, hr

if __name__ == '__main__':
    out = {'demand': dict(engines=round(D_ENG, 2), dual=round(D_DUAL, 2), dual_45=round(dem_dual(E_ANCHOR, .45), 2),
                          dual_65=round(dem_dual(E_ANCHOR, .65), 2)),
           'generic': {}, 'named': {}, 'thresholds': {}, 'sens': {}, 'buyer': {}, 'critic': {}}
    CASES = {'anchor_eng': (D_ENG, 1), 'anchor_dual': (D_DUAL, 1), 'cluster4': (D_ENG+3*0.5, 4), 'cluster_big': (6.0, 8)}
    for Q in (1, 5, 15):
        for km in (20, 26):
            for cn, (dm, nb) in CASES.items():
                dm_eff = dm if Q > 1 else min(dm, Q*0.78*0.92*0.945)      # на 1 MMscf/д покупатель ~3 МВт [О]
                for sc in ('opt', 'base', 'pess'):
                    for con in (True, False):
                        out['generic'][f'Q{Q}|{km}km|{cn}|{sc}|{"GSA" if con else "pre"}'] = run(Q, dm_eff, nb, km, sc, contract=con)
    # именованные факелы: (VIIRS 2025, минимум 2023-25, км по прямой до промузла) [О: wb/nigeria_flares.csv, funnel.py]
    NAMED = {
      'Ughelli East (NEPL) -> Warri':        (6.7, 6.3, 20.8),
      'Ughelli East (NEPL) -> Ughelli':      (6.7, 6.3, 6.5),
      'Afiesere (Shoreline) -> Ughelli':     (1.5, 1.5, 7.1),
      'Eriemu (Shoreline) -> Ughelli':       (2.1, 1.1, 8.4),
      'Utorogu (NEPL) -> Ughelli':           (3.1, 3.1, 13.2),
      'Sapele (Seplat) -> Sapele':           (6.0, 6.0, 10.8),
      'Oredo (NEPL) -> Sapele':              (8.5, 8.5, 21.8),
      'Oben (Seplat) -> Sapele':             (8.7, 8.7, 24.7),
    }
    for name, (q25, qmin, km_s) in NAMED.items():
        for cn, (dm, nb) in CASES.items():
            if cn == 'cluster_big': continue
            for sc in ('opt', 'base', 'pess'):
                km_r = km_s*SC[sc]['route']
                dm_eff = min(dm, qmin*SC[sc]['lf']*SC[sc]['up']*(1-SC[sc]['shrink']))
                out['named'][f'{name}|{cn}|{sc}'] = run(qmin, dm_eff, nb, km_r, sc, contract=True)
    # пороги: цена $/MMBtu для NPV15=0 / NPV20=0; спрос MMscf/д при $6; $/км трубы при $6
    TH = [('Q5_20km_eng', 5, D_ENG, 1, 20), ('Q5_26km_eng', 5, D_ENG, 1, 26), ('Q5_26km_dual', 5, D_DUAL, 1, 26),
          ('Q5_26km_cluster4', 5, D_ENG+1.5, 4, 26), ('Q15_26km_cluster_big', 15, 6.0, 8, 26),
          ('Q1_20km', 1, 1*0.78*0.92*0.945, 1, 20), ('Q1_8km_Afiesere', 1.5, 1.5*0.78*0.92*0.945, 1, 7.1*1.3),
          ('UghelliE_Ughelli_8km_eng', 6.3, D_ENG, 1, 6.5*1.3), ('UghelliE_Warri_27km_eng', 6.3, D_ENG, 1, 20.8*1.3),
          ('Oredo_Sapele_28km_eng', 8.5, D_ENG, 1, 21.8*1.3), ('Oredo_Sapele_28km_cluster4', 8.5, D_ENG+1.5, 4, 21.8*1.3)]
    for label, Q, dm, nb, km in TH:
        for sc in ('opt', 'base', 'pess'):
            for edti in (False, True):
                k = f'{label}|{sc}|{"EDTI" if edti else "noEDTI"}'
                out['thresholds'][k] = dict(
                    p15=solve(lambda x: run(Q, dm, nb, km, sc, ov=dict(price=x), edti=edti)['npv15'], 0.5, 80),
                    p20=solve(lambda x: run(Q, dm, nb, km, sc, ov=dict(price=x), edti=edti)['npv20'], 0.5, 80),
                    dem20_at_6=solve(lambda x: run(max(Q, x*2.5), x, nb, km, sc, ov=dict(price=6.0), edti=edti)['npv20'], 0.1, 40),   # газ не ограничивает
                )
    # чувствительности: 5 MMscf/д, 26 км, якорь на двигателях, база, GSA
    base_args = (5, D_ENG, 1, 26)
    for nm, ov, kw in (('base', {}, {}), ('EDTI', {}, dict(edti=True)), ('price_8', dict(price=8.0), {}),
                       ('price_10.75', dict(price=10.75), {}), ('price_4.7_greenville', dict(price=4.7), {}),
                       ('naira_price_12pct', {}, dict(naira_price=True)), ('build_2y', dict(build=2, ramp=0.5), {}),
                       ('conv_lag_2y', {}, dict(conv_lag=2)), ('pipe_x0.6', dict(ng=0.6), {}), ('pipe_x1.5', dict(ng=1.5), {}),
                       ('cm_1.0', dict(cm=1.0), {}), ('coll_0.75', dict(coll=0.75), {}), ('decl_15', dict(decl=0.15), {}),
                       ('opex_light_local', dict(staff=0.25, sec=0.20, surv=0.005, comm=0.08, ovh=0.10, scr=1.7), {}),
                       ('gas_in_1.5', dict(gas_in=1.5*HHV), {}), ('fee_from_operator_1.0', dict(gas_in=-1.0), {})):
        r = run(*base_args, 'base', ov=ov, **kw)
        out['sens'][nm] = dict(capex=r['design']['total'], fixed=r['fixed'], ebitda_y1=r['y1']['ebitda'],
                               ebitda_y2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
    # экономика покупателя
    for pr in (4.7, 6.0, 8.0, 10.75):
        l, hr = buyer(pr)
        out['buyer'][f'engines_lcoe_at_{pr}'] = dict(usd_kwh=round(l, 4), ngn_kwh=round(l*FX), hr_mmbtu_kwh=round(hr, 5))
    # цена газа, при которой своя генерация покупателя = нашей энергии (₦350/₦411) -- «потолок со стороны покупателя»
    for T in (350, 383, 411):
        l0, hr = buyer(0.0)
        out['buyer'][f'gas_parity_vs_power_{T}'] = round((T/FX - l0)/hr, 2)
    # dual-fuel: экономия покупателя в год и окупаемость комплекта
    diesel_mmbtu = (1780/FX/0.0365, 2020/FX/0.0365)    # ₦1 780-2 020/л [В, 10.2026], 0.0365 MMBtu/л HHV [Д]
    out['buyer']['diesel_usd_mmbtu'] = [round(x, 1) for x in diesel_mmbtu]
    for subst in (0.45, 0.55, 0.65):
        g = dem_dual(E_ANCHOR, subst)*1000*365*HHV                     # MMBtu газа в год
        for pr in (6.0, 8.0):
            sav = [g/1.03*dm - g*pr for dm in diesel_mmbtu]            # диз. энергия замещена = газ/1.03
            kit = [11500*k/1e6 for k in (30*1.5, 160*2.4)]             # $30-160/кВт [В] x монтаж 1.5-2.4 [Д], 11.5 МВт
            out['buyer'][f'dual_s{subst}_p{pr}'] = dict(gas_mmbtu=round(g), saving_musd=[round(s/1e6, 2) for s in sav],
                                                        kit_musd=[round(k, 2) for k in kit],
                                                        payback_months=[round(kit[0]/(sav[1]/1e6)*12, 1), round(kit[1]/(sav[0]/1e6)*12, 1)])
    # конверсия на газовые двигатели: капитал покупателя
    out['buyer']['engines_capex_musd'] = [round(11.5*x/1e3, 1) for x in (1500, 1667)]
    # оценка критика: 1.59 MMscf/д, маржа $5.25-7.25/Mscf, без опекса
    out['critic']['margin_musd'] = [round(1.59*365*m/1e3, 2) for m in (5.25, 7.25)]
    # сравнение с «энергией якорю»: NPV газа при цене паритета покупателя с нашей энергией ₦350/₦383 и при потолке CNG
    out['compare'] = {}
    LEAN = dict(price=8.0, staff=0.25, sec=0.20, surv=0.005, comm=0.08, ovh=0.10, scr=1.7, ng=0.6, cm=2.4)
    for label, Q, dm, nb, km in (('UghelliE->Ughelli 8km', 6.3, D_ENG, 1, 6.5*1.3), ('generic 26km', 5, D_ENG, 1, 26),
                                 ('Oredo->Sapele 28km', 8.5, D_ENG, 1, 21.8*1.3), ('cluster_big Q15 26km', 15, 6.0, 8, 26)):
        row = {}
        for pr in (6.0, 8.0, 10.75, out['buyer']['gas_parity_vs_power_350'], out['buyer']['gas_parity_vs_power_383']):
            r = run(Q, dm, nb, km, 'base', ov=dict(price=pr))
            row[f'p{pr}'] = dict(capex=r['design']['total'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
        r = run(Q, dm, nb, km, 'base', ov=LEAN, edti=True)
        row['lean_p8_EDTI'] = dict(capex=r['design']['total'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
        row['lean_p20'] = solve(lambda x: run(Q, dm, nb, km, 'base', ov=dict(LEAN, price=x), edti=True)['npv20'], 0.5, 80)
        out['compare'][label] = row
    # цена для покупателя в ₦/кВт·ч на его двигателях при нашей безубыточной цене газа (база, без EDTI)
    for k in ('UghelliE_Ughelli_8km_eng|base|noEDTI', 'Q5_26km_eng|base|noEDTI', 'Q15_26km_cluster_big|base|noEDTI',
              'UghelliE_Ughelli_8km_eng|opt|noEDTI', 'Q5_26km_eng|opt|noEDTI'):
        p20 = out['thresholds'][k]['p20']
        out['compare'][f'buyer_ngn_kwh_at_BE20 {k}'] = round(buyer(p20)[0]*FX)
    json.dump(out, open(os.path.join(HERE, 'g2g_out.json'), 'w'), indent=1, ensure_ascii=False)

    def line(k, r):
        d = r['design']
        return (f"{k:52s} cap {d['total']:5.1f} (eq {d['eq']:4.1f} pipe {d['pipe']:4.1f} {d['D']}\" ${d['pipe_per_km']:.2f}M/км hp {d['hp']}) "
                f"F {r['fixed']:.2f}(x{r['opex_vs_screen']}) S2 {r['y2']['sales']:.2f} R2 {r['y2']['rev']:.2f} "
                f"E1 {r['y1']['ebitda']:5.2f} E2 {r['y2']['ebitda']:5.2f} NPV15 {r['npv15']:6.1f} NPV20 {r['npv20']:6.1f} IRR {r['irr']} PB {r['payback']}")
    L = [f"demand: {out['demand']}"]
    for sec in ('generic', 'named'):
        L.append('== '+sec)
        for k, r in out[sec].items(): L.append(line(k, r))
    L.append('== thresholds (цена $/MMBtu NPV15=0 / NPV20=0; спрос MMscf/д для NPV20=0 при $6)')
    for k, v in out['thresholds'].items(): L.append(f'{k:48s} {v}')
    L.append('== sens (Q5, 26 км, якорь на двигателях, база, GSA)')
    for k, v in out['sens'].items(): L.append(f'{k:28s} {v}')
    L.append('== buyer'); [L.append(f'{k}: {v}') for k, v in out['buyer'].items()]
    L.append(f"== critic {out['critic']}")
    L.append('== compare'); [L.append(f'{k}: {v}') for k, v in out['compare'].items()]
    txt = '\n'.join(L); print(txt)
    open(os.path.join(HERE, 'g2g_stdout.txt'), 'w').write(txt)
