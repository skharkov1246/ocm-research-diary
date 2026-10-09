# Раунд 2. Модель сделки «энергия якорному покупателю» (одна станция, один покупатель по прямой линии 33 кВ).
# Самодостаточная копия логики ../../deal/multi_model.py (site, cfg 'P', Q=5) + ../../deal/dealcore.py (run),
# с исправлениями критика (critic.md, раздел (в)). Регрессия: при OLD-параметрах воспроизводит judge_check.py.
# Метки: [П] первичный, [В] вторичный, [О] собственный расчёт, [Д] допущение (диапазон в комментарии).
import math

TAX = 0.34            # CIT 30% + development levy 4% (NTA 2025) [В]
FX = 1330.0           # ₦/$ на 10.2026, NFEM ₦1 327–1 333 [В]; тариф долларовый, курс нужен только для пересчёта ₦<->$

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.95, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

# ---------- станция: параметры (база раунда 1, сценарий 'base' multi_model.py) ----------
OLD = dict(
    Q=5.0,            # MMscf/д заявленного факела (типовая площадка)
    avail=0.85,       # доступно/заявлено в год 1 [Д 0.6-0.85, урок NGFCP]
    decl=0.06,        # падение дебита в год [Д]
    mw_per=4.64,      # МВт брутто на 1 MMscf/д -- ОШИБКА раунда 1: HHV при КПД по LHV
    aux=0.04,         # собственные нужды [Д]
    anchor=11.0,      # МВт нетто, продаваемые якорю [Д]
    LF=0.72, av=0.92, # нагрузка якоря и готовность станции [Д]
    loss=0.03,        # потери в ЛЭП 33 кВ на 25 км [Д]
    p_dir=320.0,      # тариф, ₦/кВт·ч (в долларовом эквиваленте при FX)
    bad=0.05,         # неплатежи [Д]
    gas=0.60,         # цена газа, $/MMBtu [Д; цель GSA <=1.0 всё включено]
    opfee=0.15,       # плата оператору за подготовку/подключение, $/Mscf [Д, НЕ НАЙДЕНО]
    mmbtu_mscf=1.05,  # HHV, MMBtu/Mscf [Д 1.03-1.10]
    cm=2.6,           # поправка аудита NaCN к пакету оборудования [Д 2.4-2.8]
    kw_pkg=700.0,     # $/кВт пакет генсета (скрин) [В]
    red=1.10, kwscale=0.9,   # резерв N+1 и эффект масштаба для 5 MMscf/д (SCALE[5]) [Д]
    line_km=25.0, usd_km=90e3, subst=0.6, line_mult=2.0,   # ЛЭП 33 кВ и две подстанции (SCALE[5].line=2) [Д]
    dev=0.8,          # девелопмент [Д]
    staff_core=0.45, staff_mult=1.8, sec=0.35, comm=0.15, ovh=0.20,  # $ млн/г, персонал с экспатом [Д]
    tor=0.045, ins=0.015,    # ТОиР 4.5% [Д 3.5-5%], страхование 1.5% капитала [Д]
)
# Исправления раунда 2 к станции
NEW = dict(OLD,
    mw_per=4.40,      # [О] 1.05 MMBtu/Mscf HHV x 0.903 LHV/HHV x 293.07 кВт·ч/MMBtu x 0.38 (КПД LHV) / 24 = 4.40
    avail=0.78,       # [О/В] VIIRS: средняя загрузка установки на дебите стабильного факела 0.76-0.78 (viirs_supply.py)
    decl=0.08,        # [Д ~8%/г; VIIRS: медиана 9.2%/г за 6 лет, viirs_supply.py]
)

def station(p, r1_round=False):
    G = p['Q']*p['avail']                                   # MMscf/д доступно в год 1
    g_pow = min(G, p['anchor']/(p['mw_per']*(1-p['aux'])))  # DCQ станции = газ на полную нагрузку якоря
    mw_gross = p['mw_per']*g_pow
    cm = p['cm']
    capex = dict(
        common=(1.0+0.2*5.0)*(cm/2.6) + p['dev'],           # вход газа, осушка, учёт, площадка (как у 5 MMscf/д) [Д]
        power=mw_gross*1000*p['kw_pkg']*cm*p['red']*p['kwscale']/1e6,
        line=(p['line_km']*p['usd_km']+p['subst']*1e6)*p['line_mult']/1e6)
    C = sum(capex.values())
    if r1_round: C = round(C, 1)
    kwh_per_mscf = p['mw_per']*(1-p['aux'])*(1-p['loss'])*24     # кВт·ч у покупателя на 1 Mscf сожжённого газа
    rows = []
    for t in range(1, 11):
        Gt = G*(1-p['decl'])**(t-1)
        g_paid = min(Gt, g_pow)                              # take-or-pay на законтрактованный объём (DCQ), не выше доступного
        gas_cost = g_paid*1000*365*(p['gas']*p['mmbtu_mscf']+p['opfee'])/1e6
        gp = min(g_pow, Gt)
        dir_mw = min(p['mw_per']*gp*(1-p['aux']), p['anchor'])
        gwh = dir_mw*8760*p['av']*p['LF']*(1-p['loss'])/1e3
        rev = gwh*p['p_dir']*(1-p['bad'])/FX                 # $ млн
        fixed = (p['staff_core']+p['sec']+p['comm']+p['ovh'])*p['staff_mult'] + C*(p['tor']+p['ins'])
        if r1_round:   # как в multi_model: rev, gas, opex округлены до 0.01, капитал в fixed не округлён
            rev = round(rev, 2); fixed = round(gas_cost+fixed, 2)-round(gas_cost, 2); gas_cost = round(gas_cost, 2)
        rows.append(dict(t=t, G=Gt, g_paid=g_paid, gwh=gwh, rev=rev, gas=gas_cost, fixed=fixed))
    return dict(C=C, capex=capex, G=G, g_pow=g_pow, mw_gross=mw_gross, kwh_per_mscf=kwh_per_mscf, rows=rows)

# ---------- структуры договора ----------
# up: доля года без простоя добычи (1 мес/г простоя [Д]); gas_refund: доля оплаты газа, возвращаемая за простой;
# wc: дни дебиторки; bad: неплатежи.
STRUCTS = {
    'none':      dict(up=11/12, gas_refund=0.0, wc=60, bad=0.05),
    'esc':       dict(up=11/12, gas_refund=0.0, wc=30, bad=0.02),
    'dopgas':    dict(up=11/12, gas_refund=1.0, wc=30, bad=0.02),
    'dopmargin': dict(up=1.0,   gas_refund=0.0, wc=30, bad=0.02),
}
SLABEL = {
    'none': 'Ни DoP, ни эскроу',
    'esc': 'Только эскроу (DoP нет)',
    'dopgas': 'DoP возвращает оплату газа + эскроу (БАЗА р.2)',
    'dopmargin': 'DoP компенсирует маржу + эскроу (база судьи)',
}
# сроки и деньги: тариф в долларах (оплата в найре по NFEM), стройка 1 г, разгон 80%, инфляция долл. опекса 3%
TL_OLD = dict(dep=0.10, idx=1.0, lag=0, usd_inf=0.03, ramp=0.8, delay=1, fx_lag_q=0.0, fx_pay=False)
TL_NEW = dict(TL_OLD, fx_lag_q=0.25, fx_pay=True)   # [Д] пересчёт в найру раз в квартал + дебиторка wc дней: потеря dep*(q/2+wc/365)

def run(p, s, tl, years=10, ev=None, edti=True, depyrs=5, debt=None, tariff_cut=None, earnout=None,
        gas_extra=0.0, capex_split=None, r1_round=False, top=None):
    """p: параметры станции; s: структура; tl: сроки/валюта.
    earnout: (ставка $/Mscf реально сожжённого газа, первый операционный год) -- путь А (СП с держателем NGFCP).
    capex_split: доли капитала по годам стройки (по умолчанию весь капитал в t=0, как в раунде 1).
    top: доля DCQ под take-or-pay газа (None = 100% DCQ, как в раунде 1); оплачивается max(top*DCQ, сожжённый газ)."""
    p = dict(p, bad=s['bad'])
    st_ = station(p, r1_round); C = st_['C']; rows = st_['rows']; delay = tl['delay']
    dep = tl['dep']; wc = s['wc']
    if capex_split: cf = [-C*w for w in capex_split] + [0.0]*(delay+1-len(capex_split))
    else: cf = [-C] + [0.0]*delay
    prev = 0.0; out = []
    fxh = 1 - dep*(tl['fx_lag_q']/2 + (wc/365 if tl['fx_pay'] else 0)) if tl['idx'] >= 1 else 1.0
    for i, row in enumerate(rows):
        t = i+1; T = t+delay
        fx = (1+dep)**T; k = max(0, T-tl['lag'])
        pidx = (1+tl['idx']*((1+dep)**k-1))/fx*fxh
        if tariff_cut and t >= tariff_cut[0]: pidx *= tariff_cut[1]
        u_up = s['up']*(ev.get(t, 1.0) if ev else 1.0)
        u = u_up*(tl['ramp'] if t == 1 else 1.0)
        rev = row['rev']*pidx*u
        gas = row['gas']*(1+0.02)**(T-1)*(1-(1-s['up'])*s['gas_refund'])   # GSA с индексацией ~2%/г [Д]
        if top is not None:
            burned_d = row['gwh']*u*1e6/st_['kwh_per_mscf']/365e3
            gas *= max(top*row['g_paid'], burned_d)/row['g_paid']
        fixed = (row['fixed']+gas_extra)*(1+tl['usd_inf'])**(T-1)
        eo = 0.0
        if earnout and t >= earnout[1]:
            burned_mscf = row['gwh']*u*1e6/st_['kwh_per_mscf']
            eo = earnout[0]*burned_mscf/1e6
        e = rev-gas-fixed-eo
        d_ = C/depyrs if t <= depyrs else 0
        tx = max(0, (e-d_)*TAX); ed = (min(tx, 0.05*C) if t <= 5 else 0) if edti else 0
        w = wc/365*(rev-prev); prev = rev
        c = e-tx+ed-w+(wc/365*rev if t == years else 0)
        cf.append(c); out.append(dict(t=t, rev=rev, ebitda=e, cfads=e-tx+ed, eo=eo,
                                      burned=row['gwh']*u*1e6/st_['kwh_per_mscf']/365e3))
    r_ = irr(cf)
    res = dict(capex=round(C, 2), y1=round(out[0]['ebitda'], 2), y2=round(out[1]['ebitda'], 2),
               npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2),
               irr=None if r_ is None else round(r_*100, 1), cf=cf, rows=out, st=st_)
    if debt:
        sh, rate, ten = debt; D = C*sh; ann = D*rate/(1-(1+rate)**-ten)
        res['min_dscr'] = round(min(o['cfads'] for o in out[:ten])/ann, 2)
    return res

def be(p, s, tl, key='npv20', lo=100.0, hi=1500.0, **kw):
    """безубыточный тариф ₦/кВт·ч (долларовый эквивалент при FX)"""
    for _ in range(45):
        m = (lo+hi)/2
        if run(dict(p, p_dir=m), s, tl, **kw)[key] < 0: lo = m
        else: hi = m
    return hi

def site(name):
    """площадки матрицы"""
    if name == '25km':      # типовой факел 5 MMscf/д (покрывает Ughelli East 6.3-6.7), трасса 25 км
        return dict(NEW)
    if name == 'oredo30':   # Oredo (NEPL): VIIRS 2025 15.5, минимум 2023-25 8.5 MMscf/д [О]; трасса 30 км (28-33) [О/Д]
        return dict(NEW, Q=8.5, line_km=30.0, loss=0.03*30/25)
    raise KeyError(name)
