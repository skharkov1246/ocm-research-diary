# Скептик раунда 2 к варианту «газ генераторам покупателя» (g2g.py). Модель не меняем -- импортируем и стрессуем.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
import importlib.util, os, json, copy
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('g', os.path.join(HERE, '..', 'gas_to_genset', 'g2g.py'))
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
out = {}
U8 = (6.3, g.D_ENG, 1, 6.5*1.3); K26 = (5, g.D_ENG, 1, 26); BIG = (15, 6.0, 8, 26)
LEAN = dict(price=8.0, staff=0.25, sec=0.20, surv=0.005, comm=0.08, ovh=0.10, scr=1.7, ng=0.6, cm=2.4)

def be(args, sc='base', ov=None, key='npv20', edti=False):
    return g.solve(lambda x: g.run(*args, sc, ov=dict(ov or {}, price=x), edti=edti)[key], 0.5, 80)

# 1. Воспроизведение
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    r = g.run(*a, 'base')
    out[f'repro_{nm}'] = dict(cap=r['design']['total'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'],
                             be15=be(a, key='npv15'), be20=be(a))

# 2. Пик отбора: двигатели 11.5 МВт брутто на полной нагрузке = 11.5/4.40 = 2.61 MMscf/д товарного газа,
#    средний 1.67 -> пик/средний = 1/(0.72*0.92) = 1.51, а в модели проект = вход x1.3 (фирменная мощность 3x50% = Qdes).
#    Даём x1.51 через увеличение удельной мощности компрессии и пакетов: эквивалент -- Qdes x 1.51/1.3.
orig_design = g.design
def design_peak(Q, dem, nb, km, p, k=1.51/1.3):
    d = orig_design(Q, dem*k, nb, km, p)
    return d
# 3. Жирный газ: снижение точки росы по УВ уносит C3+ -- усушка объёма +2..6% [Д], HHV товарного газа 1.03 вместо 1.05 [Д]
def stress(args, sc='base', ov=None, peak=True, hcdp_shrink=0.03, hhv=1.03, edti=False):
    g.design = design_peak if peak else orig_design
    oldH = g.HHV; g.HHV = hhv
    p = dict(ov or {}); p['shrink'] = g.SC[sc]['shrink'] + hcdp_shrink
    r = g.run(*args, sc, ov=p, edti=edti)
    b20 = g.solve(lambda x: g.run(*args, sc, ov=dict(p, price=x), edti=edti)['npv20'], 0.5, 80)
    b15 = g.solve(lambda x: g.run(*args, sc, ov=dict(p, price=x), edti=edti)['npv15'], 0.5, 80)
    g.design = orig_design; g.HHV = oldH
    return dict(cap=r['design']['total'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], be15=b15, be20=b20)
for nm, a in (('U8', U8), ('K26', K26), ('BIG', BIG)):
    out[f'peak_hcdp_{nm}_base'] = stress(a)
    out[f'peak_hcdp_{nm}_base_hcdp6'] = stress(a, hcdp_shrink=0.06)
out['peak_hcdp_U8_lean_EDTI'] = stress(U8, ov=LEAN, edti=True)
out['peak_hcdp_U8_opt'] = stress(U8, sc='opt', hcdp_shrink=0.02, hhv=1.04)

# 4. Конкурент: отвод от магистрали (ELPS/NGMC, схема Transit Gas -> Rite Foods 18 км) тому же покупателю.
#    Газ уже подготовлен и под давлением: нет компрессии, TEG, снижения точки росы, усушки; покупка газа по
#    регулируемому ориентиру «коммерческий сектор» $2.68/MMBtu [В, Leadership 04.2026] + транспорт $0.8 (0.5-1.3) [Д].
#    Опекс тот же (персонал, охрана трассы), т.е. конкурент НЕ получает от нас скидки на нигерийские условия.
def competitor(km, gas_mmbtu, sc='base'):
    ov = dict(usd_hp=0.0, dehy=0.0, hcdp=0.0, shrink=0.0, gas_in=gas_mmbtu*g.HHV, top=0.7, p_suc=34.0)
    a = (50, g.D_ENG, 1, km)    # газ магистрали не ограничивает
    old = g.SC[sc]['decl']; g.SC[sc]['decl'] = 0.0
    r = g.run(*a, sc, ov=ov)
    b = g.solve(lambda x: g.run(*a, sc, ov=dict(ov, price=x))['npv20'], 0.5, 80)
    g.SC[sc]['decl'] = old
    return dict(cap=r['design']['total'], npv20_at6=r['npv20'], be20=b)
for km in (8.45, 18.0, 26.0):
    for gm in (2.68+0.5, 2.68+0.8, 2.68+1.3):
        out[f'competitor_spur_{km}km_gas{gm:.2f}'] = competitor(km, gm)

# 5. Системный капитал на те же 61.9 ГВт·ч/г: наш + двигатели покупателя против «энергии якорю»
eng_buyer = (11.5*1500/1e3, 11.5*1667/1e3)                  # [В] $1 500-1 667/кВт
eng_audit = 11.5*700*2.6*1.1*0.9/1e3                         # как генерация в deal_r2: $700/кВт x2.6 x1.1 x0.9
out['system_capital'] = dict(gas_8km=[round(out['repro_U8']['cap']+x, 1) for x in eng_buyer],
                             gas_26km=[round(out['repro_K26']['cap']+x, 1) for x in eng_buyer],
                             gas_26km_engines_audited=round(out['repro_K26']['cap']+eng_audit, 1),
                             power_25km=29.15, power_8km_approx=round(29.15-5.7+(8.45*90e3+0.6e6)*2/1e6, 1))

# 6. Паритет покупателя «наша энергия по ₦350/383» при согласованных затратах покупателя:
#    капитал двигателей как у нашей генерации после аудита ($1 802/кВт), ТОиР 4.5% + страховка 1.5% + персонал $0.4 млн/г [Д]
cap_kw = 700*2.6*1.1*0.9
om_kwh = (11.5*cap_kw/1e3*0.06 + 0.4)/61.9                  # $/кВт·ч
for T in (350, 383):
    l0, hr = g.buyer(0.0, capex_kw=cap_kw, om=om_kwh)
    pp = (T/g.FX - l0)/hr
    r8 = g.run(*U8, 'base', ov=dict(price=pp)); r26 = g.run(*K26, 'base', ov=dict(price=pp))
    out[f'parity_consistent_{T}'] = dict(price=round(pp, 2), om_kwh=round(om_kwh, 4),
                                         npv15_8km=r8['npv15'], npv15_26km=r26['npv15'], npv20_26km=r26['npv20'])

# 7. Потолок рынка: NPV при CNG-потолке $10.75 и при трубном газе-конкуренте
for nm, a in (('U8', U8), ('K26', K26)):
    for pr in (4.0, 4.7, 8.0, 10.75):
        r = g.run(*a, 'base', ov=dict(price=pr))
        out[f'npv_at_{pr}_{nm}'] = dict(npv15=r['npv15'], npv20=r['npv20'])
# 8. Критерий прохода «$8.2 при трубе <= $0.3 млн/км и местной команде» под стрессом пика и жирного газа
out['lean_U8_be20_EDTI_model'] = g.solve(lambda x: g.run(*U8, 'base', ov=dict(LEAN, price=x), edti=True)['npv20'], 0.5, 80)
out['lean_U8_be20_EDTI_stressed'] = out['peak_hcdp_U8_lean_EDTI']['be20']
# 9. Спрос до контракта: якорь /2.0 (база аудита) на 8 км
r = g.run(*U8, 'base', contract=False); out['U8_pre_contract'] = dict(npv15=r['npv15'], npv20=r['npv20'])
json.dump(out, open(os.path.join(HERE, 'g2g_skeptic_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)

# 10. Конкурент-инкумбент (Axxela/Transit Gas, NGMC, SNG): отвод добавляется к действующей сети --
#     предельный опекс: персонал/накладные $0.15 млн, охрана $0.10 + $6 тыс./км, община $0.05, без пола «скрин x2» [Д]
def competitor_inc(km, gas_mmbtu, cm=2.6):
    ov = dict(usd_hp=0.0, dehy=0.0, hcdp=0.0, shrink=0.0, gas_in=gas_mmbtu*g.HHV, top=0.7, p_suc=34.0,
              staff=0.10, ovh=0.05, sec=0.10, surv=0.006, comm=0.05, scr=0.0, cm=cm, dev=0.4)
    a = (50, g.D_ENG, 1, km); old = g.SC['base']['decl']; g.SC['base']['decl'] = 0.0
    r = g.run(*a, 'base', ov=ov)
    b = g.solve(lambda x: g.run(*a, 'base', ov=dict(ov, price=x))['npv20'], 0.5, 80)
    g.SC['base']['decl'] = old
    return dict(cap=r['design']['total'], fixed=r['fixed'], be20=b)
out2 = {}
for km in (8.45, 18.0, 26.0):
    for gm in (3.18, 3.98):
        out2[f'incumbent_spur_{km}km_gas{gm}'] = competitor_inc(km, gm)
# наш вариант с тем же предельным опексом (как если бы оператор NEPL/Seplat делал сам на своей площадке)
ovop = dict(staff=0.10, ovh=0.05, sec=0.10, surv=0.006, comm=0.05, scr=0.0)
out2['operator_self_U8_be20'] = g.solve(lambda x: g.run(*U8, 'base', ov=dict(ovop, price=x))['npv20'], 0.5, 80)
out2['operator_self_K26_be20'] = g.solve(lambda x: g.run(*K26, 'base', ov=dict(ovop, price=x))['npv20'], 0.5, 80)
json.dump(dict(out, **out2), open(os.path.join(HERE, 'g2g_skeptic_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out2.items(): print(k, v)
