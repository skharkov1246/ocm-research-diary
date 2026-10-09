# Пиролиз метана в расплаве (Ni-Bi) -> твёрдый углерод (товар) + H2 -> газопоршневой генератор -> э/э промпотребителю.
# Нигерия, факел 1/5/15 MMscf/д. Сопоставимо с power/captive.py (те же цены э/э, LF, линия, налог 34%, 10 лет).
# Метки входов: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение.
import json, math, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power')
FX = 1330.0
CH4_T_Y = 7300.0                 # т CH4/год на 1 MMscf/д (корпус econ_vs_mining) [О корпус]
MOL_S = CH4_T_Y * 1e6 / 16.04 / (365 * 86400)   # моль CH4/с  [О]
LHV_CH4, LHV_H2 = 802.3, 241.8   # кДж/моль [П справочник]
MW_DIRECT = 4.64                 # МВт эл. прямого генсета на 1 MMscf/д (корпус, КПД 0.38) [О корпус]
MMBTU_Y = 365 * 1000 * 1.05      # MMBtu/год на 1 MMscf/д [Д 1.03-1.1]
TAX = 0.34                       # CIT 30% + dev levy 4% (NTA 2025) [В]
YEARS = 10

def npv(r, cf): return sum(c / (1 + r) ** t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.99, 3.0
    if npv(lo + 1e-9, cf) * npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo + hi) / 2
        if npv(lo, cf) * npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

S = {
 # X — конверсия CH4; rec — извлечение углерода; heat — подвод тепла, кДж(LHV топлива)/моль CH4 (Hazer 8-10 кВт*ч/кг H2 ~ 115-145) [В]+[Д]
 'low':  dict(X=0.75, rec=0.88, wet=1.00, heat=190, eff=0.36, parasit=0.12, avail=0.75,
              c_net=100, dom_t=0,    dom_p=0,   cred_claim=0.0, cred_p=0,  slip=0.02,
              gas=1.0, pyro_capex=24.0, melt_t=80, bi=44, purif=1.0, h2_prem=1.35,
              usd_kw=1960, red=1.25, cond=1.5, line_km=40, usd_km=120e3, subst=0.8, dev=1.5,
              LF=0.60, price_ngn=250, baddebt=0.10, loss=0.04, decl=0.10, salv=0.0,
              tor=0.05, ins=0.015, staff=1.05, sec=0.5, bi_loss=2.0, cons=15, bag=25),
 'base': dict(X=0.87, rec=0.93, wet=1.00, heat=160, eff=0.38, parasit=0.10, avail=0.85,
              c_net=200, dom_t=500,  dom_p=500, cred_claim=0.5, cred_p=5,  slip=0.04,
              gas=0.5, pyro_capex=11.2, melt_t=40, bi=38, purif=0.6, h2_prem=1.25,
              usd_kw=1700, red=1.25, cond=1.0, line_km=25, usd_km=90e3, subst=0.6, dev=1.0,
              LF=0.72, price_ngn=320, baddebt=0.05, loss=0.03, decl=0.05, salv=0.15,
              tor=0.04, ins=0.012, staff=0.80, sec=0.35, bi_loss=0.5, cons=10, bag=20),
 'high': dict(X=0.95, rec=0.97, wet=1.10, heat=130, eff=0.40, parasit=0.08, avail=0.92,
              c_net=350, dom_t=1500, dom_p=800, cred_claim=1.0, cred_p=12, slip=0.06,
              gas=0.25, pyro_capex=9.6, melt_t=20, bi=20, purif=0.3, h2_prem=1.15,
              usd_kw=1500, red=1.20, cond=0.6, line_km=12, usd_km=70e3, subst=0.5, dev=0.7,
              LF=0.85, price_ngn=400, baddebt=0.03, loss=0.025, decl=0.0, salv=0.25,
              tor=0.035, ins=0.010, staff=0.60, sec=0.25, bi_loss=0.1, cons=5, bag=15),
}
SCALE = {1: dict(pyro=1.0, kw=1.0, line=1.0, staff=1.0, ds=None, disco=None),
         5: dict(pyro=5 ** 0.8 / 5, kw=0.9, line=2.0, staff=1.8),
         15: dict(pyro=15 ** 0.8 / 15, kw=0.8, line=4.0, staff=2.5)}
DS = {'low': {5: (0.3, 100), 15: (0.1, 100)}, 'base': {5: (0.5, 120), 15: (0.2, 120)}, 'high': {5: (0.7, 150), 15: (0.35, 150)}}

def balance(p):
    X = p['X']
    h2_kj = 2 * X * LHV_H2 * (1 + 0.04 * (p['wet'] - 1) / 0.1)  # жирный газ: +~4% H2 на +10% C  [О грубо]
    ch4_kj = (1 - X) * LHV_CH4
    fuel_to_gen = max(0.0, h2_kj + ch4_kj - p['heat'])          # кДж/моль CH4 на генератор
    mw_gross = MW_DIRECT * fuel_to_gen / LHV_CH4 * p['eff'] / 0.38
    mw_net = mw_gross * (1 - p['parasit'])
    c_t = CH4_T_Y * 0.75 * X * p['rec'] * p['wet']             # т C/год при 100% готовности
    h2_t = CH4_T_Y * 0.25 * X
    return dict(fuel_to_gen_kj=round(fuel_to_gen, 1), share_of_direct=round(fuel_to_gen / LHV_CH4 * p['eff'] / 0.38, 3),
                mw_gross=round(mw_gross, 2), mw_net=round(mw_net, 2), c_t_y_full=round(c_t), h2_t_y=round(h2_t),
                h2_burn_share_for_heat=round(min(1, p['heat'] / h2_kj), 2))

def run(sc, q=1, power=True, years=YEARS, over=None):
    p = dict(S[sc]); p.update(over or {})
    k = SCALE[q]; b = balance(p)
    ds, disco = (1.0, 120) if q == 1 else DS[sc][q]
    # ---- капитал
    pyro = p['pyro_capex'] * q * k['pyro']
    melt = p['melt_t'] * q * k['pyro'] * p['bi'] / 1000
    purif = p['purif'] * q * k['pyro']
    if power:
        gen = b['mw_gross'] * q * 1000 * p['usd_kw'] * k['kw'] * p['red'] * p['h2_prem'] / 1e6
        cond = p['cond'] + 0.2 * q
        line = (p['line_km'] * p['usd_km'] / 1e6 + p['subst']) * k['line']
    else:
        gen = 0.6 * q * k['kw'] * p['usd_kw'] * p['red'] / 1000  # малый генсет собственных нужд 0.6 МВт/MMscf/д [Д]
        cond = p['cond'] + 0.2 * q; line = 0.0
    capex = pyro + melt + purif + gen + cond + line + p['dev']
    # ---- поток
    price = (ds * p['price_ngn'] + (1 - ds) * disco) / FX
    LFeff = ds * p['LF'] + (1 - ds) * 0.9
    cf = [-capex]; rows = []
    for t in range(1, years + 1):
        g = (1 - p['decl']) ** (t - 1)
        ramp = 0.6 if t == 1 else 1.0                         # FOAK-разгон первого года [Д]
        a = p['avail'] * ramp * g
        c_t = b['c_t_y_full'] * q * a
        dom = min(c_t, p['dom_t'])
        rev_c = (dom * p['dom_p'] + (c_t - dom) * p['c_net']) / 1e6
        if power:
            gwh = 8760 * b['mw_net'] * q * a * min(LFeff, 1.0) * (1 - p['loss']) / 1e3
            rev_e = gwh * 1e6 * ((ds * p['price_ngn'] * (1 - p['baddebt']) + (1 - ds) * disco * 0.85) / FX) / 1e6
        else:
            gwh = 0.0; rev_e = 0.0
        # кредиты: удержанный в твёрдом виде C (доля claim) + сокращение проскока CH4 факела (GWP100=28)
        co2e = c_t * 44 / 12 * p['cred_claim'] + CH4_T_Y * q * a * p['slip'] * 28
        rev_cr = co2e * p['cred_p'] / 1e6
        mrv = 0.06 * min(1, q) if p['cred_p'] > 0 else 0.0
        gas_c = MMBTU_Y * q * a * p['gas'] / 1e6
        var = c_t * (p['bi_loss'] * p['bi'] + p['cons'] + p['bag']) / 1e6
        fixed = p['tor'] * capex + p['ins'] * capex + (p['staff'] + p['sec']) * k['staff']
        opex = gas_c + var + fixed + mrv
        rev = rev_c + rev_e + rev_cr
        e = rev - opex
        tax = max(0, TAX * (e - capex / years))
        cf.append(e - tax + (p['salv'] * capex if t == years else 0))
        rows.append(dict(t=t, c_t=round(c_t), gwh=round(gwh, 1), rev_c=round(rev_c, 2), rev_e=round(rev_e, 2),
                         rev_cr=round(rev_cr, 2), co2e_kt=round(co2e / 1e3, 1), opex=round(opex, 2), ebitda=round(e, 2)))
    y2 = rows[1]
    acc = -capex; pb = None
    for t, c in enumerate(cf[1:], 1):
        acc += c
        if acc >= 0 and pb is None: pb = t
    return dict(sc=sc, q=q, power=power, capex=round(capex, 2),
                capex_parts=dict(pyro=round(pyro, 2), melt=round(melt, 2), purif=round(purif, 2), gen=round(gen, 2),
                                 cond=round(cond, 2), line=round(line, 2), dev=p['dev']),
                balance=b, price_usd_kwh=round(price, 3), y2=y2, ebitda_y2=y2['ebitda'],
                simple_payback=round(capex / y2['ebitda'], 1) if y2['ebitda'] > 0 else None,
                dcf_payback_yr=pb, npv15=round(npv(0.15, cf), 2), npv20=round(npv(0.20, cf), 2),
                irr=(round(irr(cf) * 100, 1) if irr(cf) is not None else None))

out = {}
for sc in S:
    for q in (1, 5, 15):
        out[f'{q}MMscfd_{sc}'] = run(sc, q)
    out[f'1MMscfd_{sc}_carbon_only'] = run(sc, 1, power=False)
# чувствительности (1 MMscf/д, base)
sens = {}
for nm, o in [('c_net_0', dict(c_net=0, dom_t=0)), ('c_net_100', dict(c_net=100)), ('c_net_400_corpus', dict(c_net=400)),
              ('c_net_1000_quality', dict(c_net=1000)), ('c_net_1500_quality', dict(c_net=1500)),
              ('X_0.75', dict(X=0.75)), ('X_0.95', dict(X=0.95)), ('heat_130', dict(heat=130)), ('heat_190', dict(heat=190)),
              ('capex_pyro_24_monolith', dict(pyro_capex=24.0)), ('capex_pyro_9.6', dict(pyro_capex=9.6)),
              ('price_210_BandA', dict(price_ngn=210)), ('price_450', dict(price_ngn=450)),
              ('cred_12_full', dict(cred_claim=1.0, cred_p=12)), ('cred_30_full', dict(cred_claim=1.0, cred_p=30)), ('cred_0', dict(cred_p=0)),
              ('bi_loss_5kg_t', dict(bi_loss=5.0)), ('NO_AUDIT_pyro4.0', dict(pyro_capex=4.0, melt_t=0, purif=0, staff=0.35, sec=0.35, X=1.0, rec=1.0, c_net=400, dom_t=0)), ('NO_AUDIT_pyro4.0_carbonX087', dict(pyro_capex=4.0, melt_t=0, purif=0, staff=0.35, c_net=400, dom_t=0)), ('avail_0.7', dict(avail=0.7)), ('gas_DBP_2.18', dict(gas=2.18))]:
    r = run('base', 1, over=o)
    sens[nm] = {kk: r[kk] for kk in ['capex', 'ebitda_y2', 'simple_payback', 'npv15', 'npv20', 'irr']}
out['sens_1MMscfd_base'] = sens
# break-even цены углерода (NPV15=0) для base
def be_c(sc, q=1, power=True, r='npv15'):
    lo, hi = -500, 6000
    for _ in range(60):
        m = (lo + hi) / 2
        v = run(sc, q, power, over=dict(c_net=m, dom_t=0))[r]
        if v < 0: lo = m
        else: hi = m
    return round(m)
out['breakeven_carbon_usd_t'] = {f'{sc}_{q}': be_c(sc, q) for sc in S for q in (1, 5, 15)}
out['breakeven_carbon_usd_t_carbon_only'] = {sc: be_c(sc, 1, False) for sc in S}
cap=json.load(open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/power/captive_out.json'))
out['delta_vs_captive_npv15']={k:round(out[k]['npv15']-cap[k]['npv15'],1) for k in cap if k in out and 'npv15' in cap[k]}
out['carbon_only_cash_cost_usd_t']={sc:round(out[f'1MMscfd_{sc}_carbon_only']['y2']['opex']*1e6/out[f'1MMscfd_{sc}_carbon_only']['y2']['c_t']) for sc in S}
out['breakeven_carbon_usd_t_npv20'] = {f'{sc}_1': be_c(sc, 1, True, 'npv20') for sc in S}
# корпусный скрин без поправок и его пересчёт
scr_c = 0.75 * 7300 * 400 / 1e6
out['corpus_screen'] = dict(gross_HC=5.84, ebitda_HC=5.34, capex=4.0, carbon_only_rev=round(scr_c, 2), carbon_only_ebitda=round(scr_c - 0.5, 2),
                            audit_capex=[round(4.0 * 2.4, 1), round(4.0 * 2.8, 1)], audit_opex=[0.85, 1.2])
json.dump(out, open('/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/skeptic_pyro/dummy.json', 'w'), indent=1, ensure_ascii=False)
