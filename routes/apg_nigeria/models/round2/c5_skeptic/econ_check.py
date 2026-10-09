# Скептик раунда 2, «Конденсат C5+ (JT)». Независимая денежная модель (econ.py НЕ импортируется) [О].
# Цель: проверить базу и найти лучшее для третьей стороны сочетание, при котором вердикт «не делать» ломается.
# Поправки аудита NaCN: капитал x2.4-2.8 к пакету, опекс >= скрин x1.7-2.4, персонал с экспатом, ТОиР 3.5-5%,
# загрузка 0.78, падение 8%/г, разгон 0.7, налог 34% (с переносом убытков — мягче оценки), оборотный капитал.
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))

def model(Q, y, brent=70.0, brent_y12=None, hp=0.6, lf=0.78, decl=0.08, up=0.95, ramp=0.7, cm=2.6,
          roy=0.17, loss=0.08, tariff=3.0, diff=3.0, staff=0.40, ovh=0.15, sec=0.25, comm=0.08,
          tor=0.0425, ins=0.015, opx_mult=2.0, meoh_kg=0.0, meoh_usd_t=750, dso=75, coll=0.95,
          dev=0.5, build=1, years=10, fixed_override=None):
    s = (Q/5)**0.6
    pkg = 0.45*s + 0.30*s + 0.25 + 0.30 + 0.15*s                       # тот же пакет [Д], чтобы сравнивать
    inst = pkg*cm + dev; C = inst*(1 + 0.14*build/2)
    ssc = 1.0 if Q <= 1.5 else (1.25 if Q <= 7 else 1.5)
    fbu = (staff+ovh)*ssc + sec + comm + C*(tor+ins)
    fscr = (0.12*ssc + 0.05 + pkg*0.025)*opx_mult
    fixed = fixed_override if fixed_override is not None else max(fbu, fscr)
    cf = [-C/build]*build; prev = 0; loss_cf = 0; rows = []
    for t in range(1, years+1):
        px = brent_y12 if (brent_y12 and t <= 2) else brent
        net = (px - diff)*(1-roy)*(1-loss) - tariff
        G = Q*lf*up*hp*(1-decl)**(t-1)*(ramp if t == 1 else 1)
        bbl_d = G*y
        rev = bbl_d*365*net*coll/1e6
        opex = fixed*1.03**(t-1) + bbl_d*365*4.6*1.0/1e6 + G*365*meoh_kg/1000*meoh_usd_t/1e6
        e = rev - opex
        ti = e - C/10 + loss_cf                    # перенос убытков (мягче, чем в econ.py)
        tax = max(0, ti)*0.34; loss_cf = min(0, ti)
        wc = dso/365*(rev - prev); prev = rev
        cf.append(e - tax - wc + (dso/365*rev if t == years else 0))
        rows.append(dict(t=t, bbl_d=round(bbl_d, 1), rev=round(rev, 3), opex=round(opex, 3), ebitda=round(e, 3)))
    return dict(capex=round(C, 2), fixed=round(fixed, 3), y1=rows[0], y2=rows[1],
                npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2))

out = {}
# 1) повтор базы (средний газ AG2, без ингибитора 30 бар, 1.13 барр./MMscf) [О]
for Q in (1, 5, 15): out[f'base_Q{Q}'] = model(Q, 1.13)
# 2) лучшее для третьей стороны: газ ВД 100%, глубокий JT 30->3 с метанолом (AG6 7.38; NJTD 50->3 9.18),
#    спот Brent $95 первые 2 года, затем $70; низкий край опекса аудита (персонал 0.30, охрана 0.15, x1.7)
for Q in (15, 18, 35.8):
    for lab, y, mk in (('AG6_30-3', 7.38, 29.9), ('NJTD_50-3', 9.18, 17.7), ('ceiling_allC4_AG6', 10.2, 34.4)):
        out[f'best_Q{Q}_{lab}'] = model(Q, y, brent=70, brent_y12=95, hp=1.0, meoh_kg=mk, staff=0.30, sec=0.15,
                                        ovh=0.10, comm=0.05, tor=0.035, ins=0.01, opx_mult=1.7, cm=2.4)
        out[f'best_Q{Q}_{lab}_brent95flat'] = model(Q, y, brent=95, hp=1.0, meoh_kg=mk, staff=0.30, sec=0.15,
                                        ovh=0.10, comm=0.05, tor=0.035, ins=0.01, opx_mult=1.7, cm=2.4)
# 3) утверждение оценки «идеальный выход 100% C5+ и Brent $100 не покрывают опекс площадки»: проверка на 15 MMscf/д
for hp in (0.6, 1.0):
    out[f'claim_ideal_Brent100_Q15_hp{hp}'] = model(15, 9.84, brent=100, hp=hp)
    out[f'claim_ideal_Brent100_Q5_hp{hp}'] = model(5, 9.84, brent=100, hp=hp)
# 4) третья сторона как O&M-подрядчик оператора (персонал/охрана оператора), опекс $0.45 млн/г [Д 0.3-0.6]
for Q in (15, 18):
    out[f'om_contract_Q{Q}_AG6_30-3_base'] = model(Q, 7.38, hp=1.0, meoh_kg=29.9, fixed_override=0.45)
    out[f'om_contract_Q{Q}_AG2_noinh30_base'] = model(Q, 1.13, hp=1.0, fixed_override=0.45)
# 5) надстройка к якорю, оптимизм оценки (10.17 барр./MMscf = глубокий JT 50->3) без метанола в econ.py:
#    метанол 34.4 кг/MMscf x $600/т на 2.60 MMscf/д x 0.92 [О]
out['addon_opt_meoh_usd_m_per_y'] = round(2.60*0.92*365*34.4/1000*600/1e6, 3)
json.dump(out, open(os.path.join(HERE, 'econ_check_out.json'), 'w'), ensure_ascii=False, indent=1)
for k, v in out.items(): print(k, v)
