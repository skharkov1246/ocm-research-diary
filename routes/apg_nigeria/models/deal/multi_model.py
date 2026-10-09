# Многопродуктовая линия на одном факеле (Нигерия): NGL/LPG + энергия + CNG, с поправками аудита NaCN.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение.  Формулы модулей согласованы с
# соседними моделями фазы 2: ../power/captive.py, ../cng/cng_model.py, ../lpg/lpg_model.py.
import json, math, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'multi_out.json')
FX = 1330.0                     # ₦/$ [Д]
MW_PER_MMSCFD = 4.64            # корпус econ_vs_mining: 1000 Btu/scf, КПД 0.38 [О]
MMBTU_PER_MSCF = 1.05           # [Д 1.03-1.10]
MMBTU_PER_T_LPG = 47.3          # [В справочно]
SCM_PER_MMSCF = 28317           # [О]
MMBTU_PER_SCM = 0.03637         # [О]
CEPCI = 800/603.1               # 2018->2025 [Д]
TAX = 0.34                      # CIT 30% + development levy 4% (NTA 2025) [В]

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.95, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

SC = {
 # промоутерский скрин: пакет оборудования, без поправок
 'screen': dict(avail=1.0, decl=0.03, gas=0.25, opfee=0.0, cm=1.0, kw_pkg=700, red=1.0, LF=0.85, av=0.94,
                p_dir=400, p_disco=150, bad=0.03, line_km=12, usd_km=70e3, subst=0.5, dev=0.4,
                y_lpg=4.0, p_lpg=750, p_cond=400, ngl_share=0.0,
                p_cng=520, cng_util=0.90, cng_km=60, cng_cust=6000,
                staff_core=0.30, sec=0.15, comm=0.05, tor=0.03, ins=0.01, ovh=0.10, delay=0),
 # аудит NaCN, база: капитал x2.6 к пакету, выход по нижней границе, персонал с экспатом, ТОиР 4.5%
 'base':   dict(avail=0.85, decl=0.06, gas=0.60, opfee=0.15, cm=2.6, kw_pkg=700, red=None, LF=0.72, av=0.92,
                p_dir=320, p_disco=120, bad=0.05, line_km=25, usd_km=90e3, subst=0.6, dev=0.8,
                y_lpg=2.0, p_lpg=600, p_cond=400, ngl_share=0.25,
                p_cng=480, cng_util=0.65, cng_km=100, cng_cust=4500,
                staff_core=0.45, sec=0.35, comm=0.15, tor=0.045, ins=0.015, ovh=0.20, delay=0),
 # аудит NaCN, пессимизм: x2.8, бедный газ, дебит 60% заявленного (урок 1-го раунда NGFCP)
 'pess':   dict(avail=0.60, decl=0.10, gas=1.00, opfee=0.30, cm=2.8, kw_pkg=700, red=None, LF=0.60, av=0.88,
                p_dir=250, p_disco=100, bad=0.10, line_km=40, usd_km=120e3, subst=0.8, dev=1.0,
                y_lpg=1.0, p_lpg=450, p_cond=350, ngl_share=0.20,
                p_cng=330, cng_util=0.45, cng_km=150, cng_cust=3500,
                staff_core=0.60, sec=0.50, comm=0.30, tor=0.05, ins=0.02, ovh=0.30, delay=1),
}
# прямой спрос на э/э в радиусе линии (МВт нетто, якорные промпотребители) [Д]; остаток -> DisCo/embedded
DIRECT_MW = {1: 4.4, 5: 11.0, 15: 13.0}
# рынок CNG, доступный одной площадке (MMscf/д поданного газа) [Д]: нац. контрактный объём PCNGI ~20 MMscf/д [В]
CNG_CAP = {1: 1.0, 5: 1.5, 15: 2.0}
SCALE = {1: dict(red=1.25, kwscale=1.0, line=1.0, staff=1.0), 5: dict(red=1.10, kwscale=0.9, line=2.0, staff=1.8),
         15: dict(red=1.10, kwscale=0.8, line=4.0, staff=2.5)}

def site(Q, sc, cfg, years=10, extra_rev_per_mscf=0.0, toll=None, verbose=False, ov=None):
    """cfg: набор модулей из {'P','Pg','L','C'}; P = э/э только прямым клиентам, Pg = + остаток в DisCo"""
    p = dict(SC[sc]); p.update(ov or {}); s = SCALE[Q]
    G = Q*p['avail']                                  # MMscf/д реально доступно, год 1
    lpg_on = 'L' in cfg
    y_lpg = p['y_lpg'] if lpg_on else 0.0
    # энергия, уходящая с LPG+конденсатом, уменьшает газ для двигателей/CNG
    shrink = (y_lpg*1.2*MMBTU_PER_T_LPG)/(1000*MMBTU_PER_MSCF) if lpg_on else 0.0   # доля, на 1 MMscf
    G_res = G*(1-shrink)
    # распределение остаточного газа
    g_pow_max = G_res if 'Pg' in cfg else min(G_res, DIRECT_MW[Q]/(MW_PER_MMSCFD*0.96)) if 'P' in cfg else 0.0
    g_pow = g_pow_max
    g_cng = min(G_res-g_pow, CNG_CAP[Q]) if 'C' in cfg else 0.0
    own_mw = 0.0
    # ---------- CAPEX ----------
    cm = p['cm']; capex = {}
    red = s['red'] if p['red'] is None else p['red']
    common = (1.0+0.2*Q)*(1 if cm == 1 else cm/2.6) + p['dev']     # вход, осушка, учёт, площадка (как cond в captive.py) [Д]
    capex['common'] = common
    if g_pow > 0:
        mw = MW_PER_MMSCFD*g_pow
        capex['power'] = mw*1000*p['kw_pkg']*cm*red*s['kwscale']/1e6
        capex['line'] = (p['line_km']*p['usd_km']+p['subst']*1e6)*s['line']/1e6
    if lpg_on:
        ngl = 12*CEPCI*(Q/12)**0.6*cm                    # скрин $12 млн на 12 MMscf/д (2018), x CEPCI [Д]
        capex['ngl'] = ngl*(1-p['ngl_share'])            # фронт-енд (компрессия/осушка) общий с энергоблоком [Д 20-40%]
        capex['lpg_storage_loading'] = 0.3+0.05*Q*y_lpg  # [Д] сферы/пули, налив
    if g_cng > 0:
        st = 2.0*g_cng**0.7*cm
        sold_scm = g_cng*SCM_PER_MMSCF*0.88*p['cng_util']
        loads = sold_scm/7000
        trip_h = 2*p['cng_km']/35+3.5
        tractors = math.ceil(loads*trip_h/14*1.15)
        n_cust = math.ceil(sold_scm/p['cng_cust'])
        trailers = math.ceil((loads*0.5+loads*trip_h/24+n_cust)*1.15)
        capex['cng_station'] = st
        capex['cng_fleet'] = trailers*0.33+tractors*0.11+n_cust*0.15
        # компрессия CNG ~0.25 кВт*ч/scm [Д]: своя генерация
        own_mw = sold_scm*0.25/24/1000
        capex['own_gen'] = own_mw*1000*p['kw_pkg']*cm*1.25/1e6
    C = sum(capex.values())
    # ---------- годовой поток ----------
    cf = [-C] + ([0.0]*p['delay'])
    rows = []
    for t in range(1, years+1):
        g = (1-p['decl'])**(t-1)
        Gt = G*g
        # газ: take-or-pay на весь доступный объём + плата оператору за подготовку (НЕ НАЙДЕНО -> [Д])
        g_paid = Gt if lpg_on else min(Gt, g_pow+g_cng)       # take-or-pay на законтрактованный (перерабатываемый) объём
        gas_cost = g_paid*1000*365*(p['gas']*MMBTU_PER_MSCF+p['opfee'])/1e6
        rev = 0.0; var = 0.0; det = {}
        if lpg_on:
            lpg_t = y_lpg*Gt*365; cond_t = 0.2*y_lpg*Gt*365
            r = (lpg_t*p['p_lpg']+cond_t*p['p_cond'])/1e6
            det['lpg_t'] = round(lpg_t); det['rev_lpg'] = round(r, 2); rev += r
            var += lpg_t*25/1e6          # налив/логистика до депо $25/т [Д]
        if g_pow > 0:
            gp = min(g_pow, Gt*(1-shrink))
            mw_n = MW_PER_MMSCFD*gp*0.96
            dir_mw = min(mw_n, DIRECT_MW[Q])
            disco_mw = mw_n-dir_mw if 'Pg' in cfg else 0
            gwh_d = dir_mw*8760*p['av']*p['LF']*0.97/1e3
            gwh_g = disco_mw*8760*p['av']*0.9*0.97/1e3
            r = (gwh_d*p['p_dir']*(1-p['bad'])+gwh_g*p['p_disco']*(1-0.15))/FX
            det['GWh'] = round(gwh_d+gwh_g, 1); det['rev_power'] = round(r, 2); rev += r
            # неиспользованный газ энергоблока (LF<1) всё равно оплачен по take-or-pay выше
        if g_cng > 0:
            gc = min(g_cng, max(0, Gt*(1-shrink)-g_pow))
            sold = gc*SCM_PER_MMSCF*0.88*p['cng_util']*365
            r = sold*p['p_cng']/FX/1e6*(1-p['bad'])
            km_y = sold/7000*2*p['cng_km']
            var += km_y*0.95/1e6 + tractors*2*6000/1e6
            det['cng_scm_d'] = round(sold/365); det['rev_cng'] = round(r, 2); rev += r
        if extra_rev_per_mscf:
            r = (min(Gt, g_pow+g_cng) if not lpg_on else Gt)*1000*365*extra_rev_per_mscf/1e6; det['rev_operator_fee'] = round(r, 2); rev += r
        nmod = sum(1 for k in ('P', 'L', 'C') if any(k in c for c in cfg))
        staff = p['staff_core']*s['staff']*(1+0.35*(nmod-1))      # каждый доп. модуль +35% персонала [Д]
        fixed = staff+p['sec']*s['staff']+p['comm']*s['staff']+p['ovh']*s['staff']*(1+0.25*(nmod-1))+C*(p['tor']+p['ins'])
        opex = gas_cost+var+fixed
        e = rev-opex
        dep = C/5 if t <= 5 else 0
        tx = max(0, (e-dep)*TAX); edti = min(tx, 0.05*C) if t <= 5 else 0
        cf.append(e-tx+edti)
        rows.append(dict(t=t, rev=round(rev, 2), gas=round(gas_cost, 2), opex=round(opex, 2), ebitda=round(e, 2), **det))
    pb = None; cum = -C
    for i, c in enumerate(cf[1:], 1):
        if pb is None and c > 0 and cum+c >= 0: pb = i-1+(-cum)/c
        cum += c
    r_ = irr(cf)
    res = dict(Q=Q, sc=sc, cfg='+'.join(cfg), G_avail=round(G, 2), g_pow=round(g_pow, 2), g_cng=round(g_cng, 2),
               MW_gross=round(MW_PER_MMSCFD*g_pow, 1), capex=round(C, 1), capex_parts={k: round(v, 1) for k, v in capex.items()},
               y1=rows[0], ebitda_avg=round(sum(r['ebitda'] for r in rows)/len(rows), 2),
               npv15=round(npv(0.15, cf), 1), npv20=round(npv(0.20, cf), 1),
               irr=None if r_ is None else round(r_*100, 1), payback=None if pb is None else round(pb, 1))
    if verbose: res['rows'] = rows
    return res

out = {}
CFGS = [('P',), ('P', 'L'), ('P', 'L', 'C'), ('Pg', 'L'), ('Pg',), ('L', 'C'), ('L',)]
for Q in (1, 5, 15):
    for sc in ('screen', 'base', 'pess'):
        for cfg in CFGS:
            out[f'{Q}|{sc}|{"+".join(cfg)}'] = site(Q, sc, cfg)

# ---------- рычаги ----------
lev = {}
for Q in (5, 15):
    for fee in (0.5, 1.0, 1.75):        # плата оператора за «снятие факела», $/Mscf (штраф $3.50 [В], собираемость ~$1.7-1.9 [О])
        lev[f'{Q}|base|P+L|opfee{fee}'] = site(Q, 'base', ('P', 'L'), extra_rev_per_mscf=fee)
        lev[f'{Q}|pess|P+L|opfee{fee}'] = site(Q, 'pess', ('P', 'L'), extra_rev_per_mscf=fee)
out['levers'] = lev

# ---------- «скид как услуга»: толлинг для держателя разрешения ----------
def toll_needed(Q, sc, cfg, target=0.20, years=10):
    """Ставка $/Mscf поданного газа, при которой владелец скида получает IRR=target (после налога, EDTI).
    Держатель разрешения платит газ и плату оператору, ведёт сбыт и забирает продукт и рыночный риск;
    владелец скида несёт CAPEX, ТОиР, страховку и персонал площадки (без охраны/общины/сбыта)."""
    r = site(Q, sc, cfg); p = SC[sc]; s = SCALE[Q]; C = r['capex']
    nm = sum(1 for k in ('P', 'L', 'C') if any(k in c for c in cfg))
    fixed = p['staff_core']*s['staff']*(1+0.35*(nm-1)) + C*(p['tor']+p['ins'])
    def irr_at(fee):
        cf = [-C] + [0.0]*p['delay']
        for t in range(1, years+1):
            vol = Q*p['avail']*(1-p['decl'])**(t-1)*1000*365
            e = fee*vol/1e6 - fixed
            dep = C/5 if t <= 5 else 0
            tx = max(0, (e-dep)*TAX); ed = min(tx, 0.05*C) if t <= 5 else 0
            cf.append(e-tx+ed)
        return npv(target, cf)
    lo, hi = 0.0, 200.0
    for _ in range(100):
        m = (lo+hi)/2
        if irr_at(m) < 0: lo = m
        else: hi = m
    return round(m, 2), round(m*Q*p['avail']*1000*365/1e6, 2)
tolls = {}
for Q in (1, 5, 15):
    for sc in ('base', 'pess'):
        for cfg in [('Pg',), ('P', 'L'), ('L',)]:
            usd_mscf, musd = toll_needed(Q, sc, cfg)
            tolls[f'{Q}|{sc}|{"+".join(cfg)}'] = dict(toll_usd_per_mscf=usd_mscf, annual_need_musd=musd)
out['toll_for_IRR20'] = tolls

# ценность продукта на 1 Mscf поданного газа (выручка брутто), база [О]
p = SC['base']
kwh_per_mscf = MMBTU_PER_MSCF*293.07*0.38*0.96
out['gross_value_usd_per_mscf'] = dict(
    power_direct=round(kwh_per_mscf*p['p_dir']/FX*p['LF'], 2), power_direct_full_LF=round(kwh_per_mscf*p['p_dir']/FX, 2),
    power_disco=round(kwh_per_mscf*120/FX*0.85, 2),
    lpg_y2_600=round(2.0/1000*600+0.4/1000*400, 2), lpg_y4_750=round(4/1000*750+0.8/1000*400, 2),
    cng_480=round(1000/35.31*0.88*480/FX, 2),
    mining_hp27_7=round(MW_PER_MMSCFD*1000*24/1e3*0+265.1*27.7/1000, 2), mining_hp50=round(265.1*50/1000, 2),
    flare_penalty=3.50)
# выручка штрафов: N521.87 млрд / 203.97 млрд scf [В]
out['penalty_collected_usd_per_mscf'] = {fx: round(521.87e9/203.97e6/fx, 2) for fx in (1330, 1500, 1600)}

sens = {}
for Q in (5, 15):
    for lab, ov in [('rich_y4', dict(y_lpg=4.0)), ('rich_y5_p700', dict(y_lpg=5.0, p_lpg=700)),
                    ('rich_y4_cm1.6', dict(y_lpg=4.0, cm=1.6)), ('lean_y0.8', dict(y_lpg=0.8))]:
        for cfg in [('P',), ('P', 'L')]:
            sens[f'{Q}|base|{"+".join(cfg)}|{lab}'] = site(Q, 'base', cfg, ov=ov)
    for lab, ov in [('delay1', dict(delay=1)), ('delay2', dict(delay=2)), ('FX2000_unindexed', dict(p_dir=320*1330/2000, p_disco=120*1330/2000)),
                    ('price_BandA_210', dict(p_dir=210)), ('price_400', dict(p_dir=400)), ('LF0.5', dict(LF=0.5)), ('avail0.6', dict(avail=0.6)),
                    ('cm2.4', dict(cm=2.4)), ('cm2.8', dict(cm=2.8)), ('gas2.18_DBP', dict(gas=2.18, opfee=0.0))]:
        sens[f'{Q}|base|P|{lab}'] = site(Q, 'base', ('P',), ov=ov)
    for fee in (1.0, 1.75):
        sens[f'{Q}|base|P|opfee{fee}'] = site(Q, 'base', ('P',), extra_rev_per_mscf=fee)
        sens[f'{Q}|base|Pg|opfee{fee}'] = site(Q, 'base', ('Pg',), extra_rev_per_mscf=fee)
out['sens'] = sens
json.dump(out, open(OUT, 'w'), indent=1, ensure_ascii=False)
for k, v in out.items():
    if k in ('levers', 'sens', 'toll_for_IRR20', 'gross_value_usd_per_mscf', 'penalty_collected_usd_per_mscf'):
        print('==', k); [print(' ', kk, vv if not isinstance(vv, dict) or 'npv15' not in vv else {x: vv[x] for x in ('capex', 'npv15', 'npv20', 'irr', 'payback')}, vv['y1']['ebitda'] if isinstance(vv, dict) and 'y1' in vv else '') for kk, vv in v.items()] if isinstance(v, dict) and k != 'gross_value_usd_per_mscf' and k != 'penalty_collected_usd_per_mscf' else print(v)
    else:
        print(k, {x: v[x] for x in ('G_avail', 'g_pow', 'g_cng', 'MW_gross', 'capex', 'npv15', 'npv20', 'irr', 'payback')}, 'y1', v['y1'])
