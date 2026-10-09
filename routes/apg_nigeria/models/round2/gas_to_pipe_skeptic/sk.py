# Скептик к «газ в трубу» (g2p.py). Метки: [О] расчёт, [Д] допущение, [В] вторичный.
# Проверки: (1) транспортная надбавка к DBP при собственном отводе до ворот ТЭС (MYTO 2014: $0.80 за транспорт [В]);
# (2) двойной счёт простоя апстрима (VIIRS-годовой объём уже включает простой) -> up=1.0;
# (3) ТОиР/страховка на трубу 1%/0.5% вместо 4.25%/1.5% на весь капитал;
# (4) газ на входе по полу NGFCP $0.25 [В] вместо $0.78;
# (5) капитал пакета по ориентиру Seplat $0.22-0.36 млн на MMscf/д [В] x2.6 (аудит NaCN);
# (6) масштаб: какой стабильный факел нужен с этими поправками.
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gas_to_pipe'))
import g2p
TR = 0.80/1.05   # $/MMBtu транспортная составляющая, если её получает владелец отвода [В, 2014; текущий тариф NGIC НЕ НАЙДЕН]
def run_pipe_tor(Q, km, sc, ov, contract=True):
    # ТОиР+страховка: 4.25%+1.5% только на непрубную часть, на трубу 1%+0.5%
    p = dict(g2p.SC[sc]); p.update(ov)
    cx = g2p.capex(Q, km, p)
    share_pipe = cx['pipe']/cx['total']
    tor_eff = p['tor']*(1-share_pipe) + 0.01*share_pipe
    ins_eff = p['ins']*(1-share_pipe) + 0.005*share_pipe
    return g2p.run(Q, km, sc, ov=dict(ov, tor=tor_eff, ins=ins_eff), contract=contract)
cases = {
 'base': {},
 '+transport': dict(price=2.18+TR),
 '+up1.0': dict(up=1.0),
 '+gas0.25': dict(gas_in=0.25),
 '+all3': dict(price=2.18+TR, up=1.0, gas_in=0.25),
 '+all3+coll1_nofx_dso45 (эскроу)': dict(price=2.18+TR, up=1.0, gas_in=0.25, coll=1.0, fxloss=0.0, dso=45),
}
out = {}
for site, Q, km in (('Ughelli East', 6.7, 3.7), ('Afam Umuosi qmin', 2.9, 0.9), ('generic5_4km', 5, 4), ('generic15_4km', 15, 4)):
    for nm, ov in cases.items():
        r = run_pipe_tor(Q, km, 'base', ov)
        out[f'{site}|{nm}'] = dict(cap=r['capex']['total'], fixed=r['fixed_opex'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
    # капитал по ориентиру Seplat: пакет = 0.22..0.36 млн/MMscfd * Q (вместо снизу), x2.6
    for b in (0.22, 0.36):
        ov = dict(price=2.18+TR, up=1.0, gas_in=0.25)
        base_pkg = g2p.capex(Q, km, dict(g2p.SC['base'], **ov))['pkg']
        scale = b*Q/base_pkg
        ov2 = dict(ov, usd_hp=1250*scale, dehy=0.35*scale, meter=0.40*scale)
        r = run_pipe_tor(Q, km, 'base', ov2)
        out[f'{site}|all3+Seplat_capex_{b}'] = dict(cap=r['capex']['total'], fixed=r['fixed_opex'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
# масштаб с поправками all3 + эскроу, 4 км
f = lambda x: run_pipe_tor(x, 4, 'base', dict(price=2.18+TR, up=1.0, gas_in=0.25, coll=1.0, fxloss=0.0, dso=45))['npv20']
out['Q_threshold_all3_escrow_4km'] = g2p.solve(f, 0.5, 300)
f = lambda x: run_pipe_tor(x, 4, 'base', dict(price=2.18+TR, up=1.0, gas_in=0.25))['npv20']
out['Q_threshold_all3_4km'] = g2p.solve(f, 0.5, 300)
# валовая маржа на Mscf входа, база vs all3
for nm, pr, gi, up in (('base', 2.18, 0.78, 0.92), ('all3', 2.18+TR, 0.25, 1.0)):
    m = pr*1.05*(1-0.065)*0.85*(1-0.12*120/365) - gi
    out[f'margin_per_Mscf_in|{nm}'] = round(m, 2)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sk_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)
