# Раунд 2. Связка «энергия якорю + газ в трубу для остатка факела» (идея критика №5, critic.md (а)5).
# Одна площадка, одна проектная компания: газопоршневая станция 11 МВт нетто по прямой линии 33 кВ одному
# покупателю (модель ../deal/deal_r2.py) + компрессия, TEG-осушка, точка росы и отвод до соседней газовой ТЭС
# (модель ../gas_to_pipe/g2p.py) для всего газа, который станция не сжигает, включая «провал» нагрузки якоря.
# Сравниваем на ОДНОЙ временной шкале (стройка = delay лет без выручки, как в deal_r2; IDC в NPV не входит —
# NPV без долга):
#   power    — только энергия (воспроизводит deal_r2.run, регрессия ниже);
#   pipe     — только труба на весь факел (своя компания, полный опекс g2p);
#   separate — энергия + отдельная компания-труба на газ сверх DCQ станции (два ToP, два набора постоянных затрат,
#              раздельные налоги; провал нагрузки станции оплачен по ToP и сожжён);
#   bundle   — связка: общие вход/учёт/площадка/охрана/девелопмент, общий коридор, общий GSA (ToP-«хвост»
#              станции уходит в трубу), консолидированный налог.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение (диапазон в комментарии).
# Поправки аудита NaCN: капитал x2.4/2.6/2.8 к пакету (обе части), денежный опекс >= скрин x1.7/2.0/2.4,
# персонал с экспатом (в энергоблоке, staff_mult 1.8), ТОиР 3.5-5%, выход по нижней границе (усадка компрессии 6.5%,
# собственные нужды 5%), спрос /1.3-2.7 до подписания GSA трубы. Второй год: стройка 1-2 г, разгон, загрузка 0.78
# (VIIRS), падение 5-12%/г, девальвация 8-16%/г на дебиторке и квартальном лаге, оборотный капитал, налог 34%.
# Мощность 4.40 МВт на MMscf/д (исправлено, deal_r2.NEW).
# Запуск: python3 bundle.py  -> bundle_out.json и печать рядом.
import os, sys, json, copy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'deal'))
sys.path.insert(0, os.path.join(HERE, '..', 'gas_to_pipe'))
import deal_r2 as D
import g2p as G

TAX = 0.34
YEARS = 10
SHRINK = G.SHRINK + 0.01          # 5.5% топливо приводов компрессоров + 1% потерь [Д, как в g2p]
UP_PROD = 11/12                   # простой добычи 1 мес/г [Д, как deal_r2 STRUCTS]
npv, irr = D.npv, D.irr

# Сценарии. Труба — g2p.SC (opt/base/pess); энергоблок: cm 2.4/2.6/2.8, загрузка факела и падение те же, что у трубы.
SCN = {
  'opt':  dict(cm_pow=2.4, dep=0.08, bo=0.95, delay=1),   # bo — доля времени, когда ТЭС-покупатель берёт газ [Д 0.75-0.95]
  'base': dict(cm_pow=2.6, dep=0.12, bo=0.90, delay=1),   # dep — девальвация найры, %/г [Д 8-16%]
  'pess': dict(cm_pow=2.8, dep=0.16, bo=0.75, delay=2),   # delay — стройка, лет без выручки
}
# Что в связке общее (снижает приростной капитал и опекс трубы) [Д]
SHARE = dict(
  meter=0.15,        # вход/сепарация/учёт уже в «common» энергоблока; остаётся коммерческий узел учёта у ТЭС [Д 0.10-0.25]
  dev=0.30,          # изыскания/разрешения трубы при общем девелопменте [Д 0.2-0.5] (отдельно: 0.6/0.8/1.0)
  dehy_k=0.7,        # часть осушки уже в «common» энергоблока; TEG до спецификации трубы — своя [Д 0.6-0.9]
  corridor_k=0.9,    # общий коридор с ЛЭП (одно направление: Oredo->Sapele, Ughelli East->Ughelli) [Д 0.8-1.0]
  staff=0.20, ovh=0.05, sec=0.10, comm=0.05,   # приростной опекс, $ млн/г: компрессорщики без 2-го экспата,
                                               # патруль трассы, работа с общиной на трассе [Д 0.12-0.30 / 0.05-0.20]
)

def pipe_capex(Qd, km, pp, shared):
    """капитал трубного блока без IDC, $ млн; Qd — расчётный поток в трубу, MMscf/д"""
    hpq = G.hp_per_mmscfd(pp['p_suc'], pp['p_dis'])
    hp = hpq*Qd*1.5                                               # 3x50% (N+1) [Д]
    dehy = pp['dehy']*Qd**0.6*(SHARE['dehy_k'] if shared else 1.0)
    meter = (SHARE['meter'] if shared else pp['meter']) + 0.03*Qd
    pkg = hp*pp['usd_hp']/1e6 + dehy + meter
    L = km*pp['route']
    pipe = L*pp['pipe_km']*(1 + 0.15*(Qd > 7))*(SHARE['corridor_k'] if shared else 1.0)
    dev = SHARE['dev'] if shared else pp['dev']
    C = pkg*pp['cm'] + pipe + dev + 0.03                          # + лицензия NMDPRA ~$22-42 тыс. [В]
    screen = pkg + L*0.4 + 0.3                                    # промоутерский скрин [Д, как g2p]
    return dict(C=C, pkg=pkg, pipe=pipe, L=L, hp=hp, screen=screen)

def pipe_fixed(Q_site, cx, pp, shared):
    if shared:
        bu = SHARE['staff'] + SHARE['ovh'] + SHARE['sec'] + SHARE['comm'] + 0.01*cx['L'] + cx['C']*(pp['tor']+pp['ins'])
        scr = 0.08 + cx['screen']*0.025
    else:
        ss = G.staff_scale(Q_site)
        bu = (pp['staff']+pp['ovh'])*ss + pp['sec'] + 0.01*cx['L'] + pp['comm'] + cx['C']*(pp['tor']+pp['ins'])
        scr = 0.20*ss + 0.10 + cx['screen']*0.025
    return max(bu, scr*pp['opx_screen_mult'])

def model(Q, line_km, pipe_km, sc='base', T=405.0, edti=True, mode='bundle', struct='dopgas',
          contract=True, gas_allin=None, fee_op=0.0, price=None, years=YEARS, top_pipe=None, pov=None, shr=None):
    """Q — дебит факела (VIIRS), MMscf/д; T — тариф ₦/кВт·ч в долларовом эквиваленте (₦1 330/$);
    fee_op — плата оператора за снятие факела, $/Mscf газа, ушедшего в трубу (путь Б); price — цена газа у ТЭС, $/MMBtu."""
    global SHARE
    saved = SHARE
    if shr: SHARE = dict(SHARE, **shr)
    try:
        return _model(Q, line_km, pipe_km, sc, T, edti, mode, struct, contract, gas_allin, fee_op, price, years, top_pipe, pov)
    finally:
        SHARE = saved

def _model(Q, line_km, pipe_km, sc, T, edti, mode, struct, contract, gas_allin, fee_op, price, years, top_pipe, pov):
    s = D.STRUCTS[struct]; z = SCN[sc]
    pp = dict(G.SC[sc]);  pp.update(pov or {})
    if price is not None: pp['price'] = price
    delay = z['delay']; dep = z['dep']
    tl = dict(D.TL_NEW, delay=delay, dep=dep)
    depyrs = 5 if edti else 10
    p = dict(D.NEW, Q=Q, line_km=line_km, loss=0.03*line_km/25, avail=pp['lf'], decl=pp['decl'],
             cm=z['cm_pow'], p_dir=T, bad=s['bad'])
    if gas_allin is not None: p['gas'] = gas_allin - 0.15/1.05
    gas_mscf = p['gas']*p['mmbtu_mscf'] + p['opfee']           # $/Mscf, одна цена GSA для обеих частей
    has_pow = mode in ('power', 'bundle', 'separate')
    has_pipe = mode in ('pipe', 'bundle', 'separate')
    shared = mode == 'bundle'
    st = D.station(p)
    G_av = [r['G'] for r in st['rows']]                         # доступно, MMscf/д (Q*загрузка*падение)
    # --- энергоблок: построчно как deal_r2.run (TL_NEW, idx=1, lag=0) ---
    Cp = st['C'] if has_pow else 0.0
    fxh = 1 - dep*(tl['fx_lag_q']/2 + s['wc']/365)
    pw = []
    for i, row in enumerate(st['rows']):
        t = i+1; Tt = t+delay
        if not has_pow:
            pw.append(dict(rev=0, gas=0, fixed=0, burned=0, slack=0, gpaid=0)); continue
        u = s['up']*(tl['ramp'] if t == 1 else 1.0)
        rev = row['rev']*fxh*u
        gpaid = row['g_paid']*(1-(1-s['up'])*s['gas_refund'])    # оплаченный газ, MMscf/д (ToP на DCQ)
        gas = row['gas']*(1.02)**(Tt-1)*(1-(1-s['up'])*s['gas_refund'])
        fixed = row['fixed']*(1+tl['usd_inf'])**(Tt-1)
        burned = row['gwh']*u*1e6/st['kwh_per_mscf']/365e3
        pw.append(dict(rev=rev, gas=gas, fixed=fixed, burned=burned, slack=max(0.0, gpaid-burned), gpaid=gpaid))
    # --- трубный блок ---
    dem = 1.0 if contract else pp['demand']
    own = pp['own']
    if mode == 'pipe':
        base_flow = [G_av[i]*UP_PROD*own for i in range(years)]
    elif mode == 'bundle':                                       # весь газ, не сожжённый станцией
        base_flow = [max(0.0, G_av[i]*UP_PROD - pw[i]['burned'])*own for i in range(years)]
    else:                                                        # separate: только газ сверх DCQ станции
        base_flow = [max(0.0, G_av[i]-st['g_pow'])*UP_PROD*own for i in range(years)]
    Qd = base_flow[1] if has_pipe else 0.0                       # расчёт на 2-й год (после разгона станции) [Д]
    if has_pipe and Qd < 0.05: has_pipe = False
    cx = pipe_capex(Qd, pipe_km, pp, shared) if has_pipe else dict(C=0.0, L=0, screen=0, pkg=0, pipe=0, hp=0)
    Cg = cx['C']
    Fg = pipe_fixed(Q, cx, pp, shared) if has_pipe else 0.0
    cap = Qd*1.25                                                # перегруз сверх расчёта при простое станции [Д]
    cf_p = [-Cp] + [0.0]*delay; cf_g = [-Cg] + [0.0]*delay
    prev_p = prev_g = 0.0; rows = []
    for i in range(years):
        t = i+1; Tt = t+delay; a = pw[i]
        # труба
        if has_pipe:
            avail = base_flow[i]
            take = min(cap, avail*z['bo']*dem)*(pp['ramp'] if t == 1 else 1.0)
            sales = take*(1-SHRINK)
            pr = (pp['price']*(1+pp['pdrift'])**(t-1) - pp['transport'])*pp['hv']           # $/Mscf
            rev_g = sales*365e3*pr/1e6*pp['coll']*(1 - dep*pp['dso']/365)
            free = a['slack'] if shared else 0.0                 # провал нагрузки станции уже оплачен по ToP
            payable = (avail if (top_pipe if top_pipe is not None else pp['top_all']) else take)
            gas_g = max(0.0, payable - free)*365e3*gas_mscf/1e6*(1.02)**(Tt-1)
            fee = take*365e3*fee_op/1e6
            fix_g = Fg*(1.03)**(Tt-1)
            e_g = rev_g + fee - gas_g - fix_g
        else:
            take = sales = rev_g = gas_g = fee = fix_g = e_g = 0.0
        e_p = a['rev'] - a['gas'] - a['fixed']
        dp = Cp/depyrs if t <= depyrs else 0.0
        dg = Cg/depyrs if t <= depyrs else 0.0
        w_p = s['wc']/365*(a['rev']-prev_p); prev_p = a['rev']
        w_g = pp['dso']/365*(rev_g+fee-prev_g); prev_g = rev_g+fee
        end_p = s['wc']/365*a['rev'] if t == years else 0.0
        end_g = pp['dso']/365*(rev_g+fee) if t == years else 0.0
        if shared:   # консолидированный налог и EDTI на весь капитал
            tx = max(0.0, (e_p+e_g-dp-dg)*TAX)
            ed = (min(tx, 0.05*(Cp+Cg)) if t <= 5 else 0.0) if edti else 0.0
            cf_p.append(e_p+e_g - tx + ed - w_p - w_g + end_p + end_g); cf_g.append(0.0)
        else:
            tp = max(0.0, (e_p-dp)*TAX); ep_ = (min(tp, 0.05*Cp) if t <= 5 else 0.0) if edti else 0.0
            tg = max(0.0, (e_g-dg)*TAX); eg_ = (min(tg, 0.05*Cg) if t <= 5 else 0.0) if edti else 0.0
            cf_p.append(e_p - tp + ep_ - w_p + end_p); cf_g.append(e_g - tg + eg_ - w_g + end_g)
        rows.append(dict(t=t, avail=round(G_av[i]*UP_PROD, 2), burned=round(a['burned'], 2), slack=round(a['slack'], 2),
                         pipe_take=round(take, 2), pipe_sales=round(sales, 2), ebitda_pow=round(e_p, 2),
                         ebitda_pipe=round(e_g, 2), ebitda=round(e_p+e_g, 2)))
    cf = [x+y for x, y in zip(cf_p, cf_g)]
    r = irr(cf)
    acc = 0.0; pb = None
    for k, c in enumerate(cf):
        if pb is None and k > 0 and acc < 0 and acc+c >= 0: pb = round(k-1+(-acc)/c, 1)
        acc += c
    flared = sum(max(0.0, rw['avail']*own - rw['burned'] - rw['pipe_take']) if mode != 'pipe' else rw['avail']*own - rw['pipe_take']
                 for rw in rows)/years
    used = sum(rw['burned'] + rw['pipe_take'] for rw in rows)/years
    return dict(mode=mode, Q=Q, sc=sc, T=T, edti=edti, capex=round(Cp+Cg, 2), capex_pow=round(Cp, 2), capex_pipe=round(Cg, 2),
                Qd_pipe=round(Qd, 2), hp=round(cx['hp']), pipe_len=round(cx['L'], 1), fixed_pipe=round(Fg, 2),
                y1=rows[0]['ebitda'], y2=rows[1]['ebitda'], y2_pipe=rows[1]['ebitda_pipe'],
                npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2), irr=None if r is None else round(r*100, 1),
                payback_from_t0=pb, use_share=round(used/max(1e-9, used+flared), 2),
                pipe_sales_mmscfd=[rw['pipe_sales'] for rw in rows], rows=rows, cf=[round(c, 3) for c in cf])

def solve(fn, lo, hi, it=60):
    flo, fhi = fn(lo), fn(hi)
    if (flo >= 0) == (fhi >= 0): return None
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm >= 0) == (flo >= 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 2)

# ---------------- углерод AM0009 на газ, поставленный в трубу (наложение, ожидаемое значение) ----------------
# Параметры — carbon_sk3/y2.py и phase3_context.md (скептики 2-3 по углероду):
V0 = 17.7                                   # тыс. т CO2e/г на 1 MMscf/д реально поставленного газа, нижняя граница [О корпуса]
NETP = {'VCM4': 2.32, 'C10': 5.00, 'I20': 13.60}   # $/т нетто после вычетов [О корпуса]
P_OK = 0.6                                  # вероятность письма NCCC о допустимости [Д корпуса]
def carbon(sales, price, best=True, window=8, delay_build=1, r=0.15):
    dem, mrvm, capm, infl, lag = (1.3, 1.7, 2.4, 0.03, 2) if best else (2.7, 2.4, 2.8, 0.06, 3)
    cap = 0.30*capm*1.5; mrv0 = 0.08*mrvm*1.3
    n = len(sales)
    cash = [0.0]*(delay_build+n+lag+2); cash[0] = -cap
    for t in range(1, min(window, n)+1):
        T = t+delay_build
        cash[T] -= mrv0*(1+infl)**t*(1-TAX)
        cash[T+lag-1] += V0*sales[t-1]*NETP[price]/dem/1000*(1-TAX)
    ok = npv(.15, cash); ok20 = npv(.20, cash)
    ev15 = P_OK*ok - (1-P_OK)*0.5*cap; ev20 = P_OK*ok20 - (1-P_OK)*0.5*cap
    return dict(npv15_if_ok=round(ok, 2), ev15=round(ev15, 2), ev20=round(ev20, 2), kt_y1=round(V0*sales[1], 1))

SITES = {   # Q (VIIRS, MMscf/д) [О: wb/nigeria_flares.csv], ЛЭП км [Д/О судьи], отвод до ТЭС км по прямой [О g2p, GEM wiki [В]]
  'Ughelli East 2025 (6.7) -> Transcorp Ughelli':   (6.7, 25.0, 3.7),
  'Ughelli East min23-25 (6.3)':                    (6.3, 25.0, 3.7),
  'Ughelli East 6.7, ЛЭП 10 км [Д]':                 (6.7, 10.0, 3.7),
  'Oredo min23-25 (8.5) -> Sapele PS':              (8.5, 30.0, 16.4),
  'Oredo 2025 (15.5) -> Sapele PS':                 (15.5, 30.0, 16.4),
  'типовой 5 MMscf/д, отвод 4 км':                   (5.0, 25.0, 4.0),
  'типовой 15 MMscf/д, отвод 8 км':                  (15.0, 25.0, 8.0),
}

if __name__ == '__main__':
    OUT = {'regression': {}, 'sites': {}, 'be_tariff': {}, 'thresholds': {}, 'carbon': {}, 'sens': {}, 'q1': {}}
    # 0. регрессия: mode='power' воспроизводит deal_r2.run (25 км, 5 MMscf/д, база р.2)
    for T, e in ((350, True), (411, False)):
        ref = D.run(dict(D.NEW, p_dir=T), D.STRUCTS['dopgas'], dict(D.TL_NEW, dep=0.12), edti=e, depyrs=5 if e else 10)
        mine = model(5.0, 25.0, 4.0, 'base', T, e, 'power')
        OUT['regression'][f'{T}|{e}'] = dict(deal_r2=ref['npv15'], bundle_py=mine['npv15'])
    print('регрессия NPV15 deal_r2 vs bundle.py (power, base):', OUT['regression'])

    TARIFFS = (350, 377, 411, 440)
    for nm, (Q, lk, pk) in SITES.items():
        for sc in ('opt', 'base', 'pess'):
            for e in (True, False):
                for T in TARIFFS:
                    k = f'{nm}|{sc}|{"EDTI" if e else "noEDTI"}|{T}'
                    res = {m: model(Q, lk, pk, sc, T, e, m) for m in ('power', 'pipe', 'separate', 'bundle')}
                    pre = model(Q, lk, pk, sc, T, e, 'bundle', contract=False)
                    OUT['sites'][k] = {m: {kk: v[kk] for kk in ('capex', 'capex_pow', 'capex_pipe', 'Qd_pipe', 'fixed_pipe', 'y1', 'y2',
                                                               'y2_pipe', 'npv15', 'npv20', 'irr', 'payback_from_t0', 'use_share')}
                                       for m, v in res.items()}
                    OUT['sites'][k]['bundle_preGSA'] = {kk: pre[kk] for kk in ('npv15', 'npv20', 'y2')}
                    OUT['sites'][k]['increment_npv15'] = round(res['bundle']['npv15']-res['power']['npv15'], 2)
                    OUT['sites'][k]['increment_npv20'] = round(res['bundle']['npv20']-res['power']['npv20'], 2)
                    OUT['sites'][k]['synergy_vs_separate_npv15'] = round(res['bundle']['npv15']-res['separate']['npv15'], 2)
            # безубыточный тариф power vs bundle
            for e in (True, False):
                for key in ('npv15', 'npv20'):
                    for m in ('power', 'bundle'):
                        OUT['be_tariff'][f'{nm}|{sc}|{"EDTI" if e else "noEDTI"}|{key}|{m}'] = solve(
                            lambda T: model(Q, lk, pk, sc, T, e, m)[key], 150, 1200)
    # пороги для приращения «труба в связке» (NPV20 приращения = 0) при T=405 (порог LOI р.2, EDTI), база
    for nm, (Q, lk, pk) in SITES.items():
        for sc in ('opt', 'base'):
            T = 405; e = True
            p0 = model(Q, lk, pk, sc, T, e, 'power')['npv20']
            inc = lambda **kw: model(Q, lk, pk, sc, T, e, 'bundle', **kw)['npv20'] - p0
            OUT['thresholds'][f'{nm}|{sc}'] = dict(
                inc_npv20=round(inc(), 2),
                gas_sale_price_usd_mmbtu=solve(lambda x: inc(price=x), 0.5, 30),
                operator_fee_usd_mscf=solve(lambda x: inc(fee_op=x), 0.0, 15),
                gsa_allin_usd_mmbtu=solve(lambda x: inc(gas_allin=x), 0.0, 3.0),
                collection=solve(lambda x: inc(pov=dict(coll=x)), 0.3, 1.0),
                pipe_cost_usd_m_per_km=solve(lambda x: inc(pov=dict(pipe_km=x)), 0.0, 3.0))
    # углерод на трубный газ связки (EV, после налога), база, T=405 EDTI
    for nm, (Q, lk, pk) in SITES.items():
        b = model(Q, lk, pk, 'base', 405, True, 'bundle')
        o = model(Q, lk, pk, 'opt', 405, True, 'bundle')
        for pr in NETP:
            for best in (True, False):
                for w in (4, 8):
                    OUT['carbon'][f'{nm}|{pr}|{"best" if best else "worst"}|{w}y'] = dict(
                        base_flow=carbon(b['pipe_sales_mmscfd'], pr, best, w), opt_flow=carbon(o['pipe_sales_mmscfd'], pr, best, w))
    # чувствительность, Oredo 8.5 и Ughelli East 6.7, база, T=405 EDTI
    for nm in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo min23-25 (8.5) -> Sapele PS', 'Oredo 2025 (15.5) -> Sapele PS'):
        Q, lk, pk = SITES[nm]
        p0 = model(Q, lk, pk, 'base', 405, True, 'power')
        cases = {
          'база': {}, 'коммерческая цена $2.68': dict(price=2.68), 'газ NGFCP пол $0.25+0.15/Mscf': dict(gas_allin=0.40/1.05),
          'эскроу у ТЭС: сбор 0.98, 45 дн': dict(pov=dict(coll=0.98, dso=45)),
          'плата оператора $1.0/Mscf': dict(fee_op=1.0), 'плата оператора $1.75/Mscf': dict(fee_op=1.75),
          'капитал трубы x2.4 вместо x2.6': dict(pov=dict(cm=2.4)), 'труба $0.5 млн/км': dict(pov=dict(pipe_km=0.5)),
          'нет синергии (SHARE как отдельно)': dict(shr=dict(meter=0.40, dev=0.8, dehy_k=1.0, corridor_k=1.0, staff=0.55*1.3, ovh=0.2*1.3, sec=0.35, comm=0.15)),
          'сильная синергия (staff 0.12, sec 0, dev 0.2)': dict(shr=dict(staff=0.12, sec=0.0, dev=0.2, corridor_k=0.8, dehy_k=0.6)),
          'ToP трубы на весь доступный газ': dict(top_pipe=True),
          'всё вместе: $2.68 + эскроу + пол NGFCP': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45)),
        }
        for cn, kw in cases.items():
            r = model(Q, lk, pk, 'base', 405, True, 'bundle', **kw)
            OUT['sens'][f'{nm}|{cn}'] = dict(capex=r['capex'], y2=r['y2'], npv15=r['npv15'], npv20=r['npv20'],
                                             inc15=round(r['npv15']-p0['npv15'], 2), inc20=round(r['npv20']-p0['npv20'], 2), irr=r['irr'])
    # 1 MMscf/д: якоря 11 МВт нет (доступно ~0.78 MMscf/д = 3.3 МВт) -> только труба
    for sc in ('opt', 'base', 'pess'):
        for km in (1.0, 4.0):
            r = model(1.0, 25.0, km, sc, 405, True, 'pipe')
            OUT['q1'][f'{sc}|{km}km'] = dict(capex=r['capex'], y2=r['y2'], npv15=r['npv15'], npv20=r['npv20'])
    json.dump(OUT, open(os.path.join(HERE, 'bundle_out.json'), 'w'), indent=1, ensure_ascii=False)

    def show(k):
        v = OUT['sites'][k]
        print(f'  {k}')
        for m in ('power', 'pipe', 'separate', 'bundle'):
            x = v[m]
            print(f'     {m:9s} cap {x["capex"]:6.1f} (pow {x["capex_pow"]:5.1f} pipe {x["capex_pipe"]:5.1f}, Qd {x["Qd_pipe"]:5.2f}) '
                  f'E1 {x["y1"]:6.2f} E2 {x["y2"]:6.2f} (pipe {x["y2_pipe"]:5.2f}) NPV15 {x["npv15"]:6.1f} NPV20 {x["npv20"]:6.1f} '
                  f'IRR {x["irr"]} PB {x["payback_from_t0"]} газ исп. {x["use_share"]}')
        print(f'     до GSA трубы: {v["bundle_preGSA"]} | приращение трубы NPV15 {v["increment_npv15"]} NPV20 {v["increment_npv20"]}'
              f' | связка против раздельно NPV15 {v["synergy_vs_separate_npv15"]}')
    print('\n=== площадки (база р.2: DoP-газ + эскроу) ===')
    for k in OUT['sites']:
        if ('|405' in k) or ('|411' in k and '|EDTI' in k) or ('|377' in k and 'base' in k): show(k)
    print('\n=== безубыточный тариф ₦/кВт·ч: power -> bundle ===')
    for nm in SITES:
        for sc in ('opt', 'base', 'pess'):
            for e in ('EDTI', 'noEDTI'):
                a = [OUT['be_tariff'][f'{nm}|{sc}|{e}|{k}|{m}'] for k in ('npv15', 'npv20') for m in ('power', 'bundle')]
                print(f'  {nm:46s} {sc:4s} {e:6s} NPV15=0: {a[0]} -> {a[1]} | NPV20=0: {a[2]} -> {a[3]}')
    print('\n=== пороги приращения трубы (NPV20 приращения = 0), T=₦405, EDTI ===')
    for k, v in OUT['thresholds'].items(): print(' ', k, v)
    print('\n=== углерод AM0009 на газ трубы связки (база потока; EV после налога, P=0.6) ===')
    for k, v in OUT['carbon'].items():
        if '|8y' in k or 'I20' in k: print(' ', k, v)
    print('\n=== чувствительность (база, T=₦405, EDTI) ===')
    for k, v in OUT['sens'].items(): print(' ', k, v)
    print('\n=== 1 MMscf/д: только труба ===', OUT['q1'])
