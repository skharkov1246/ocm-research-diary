# Скептик BOOT: проверка ключевых чисел boot.py и вердикта «не делать».
# Все функции модели — из ../compression_boot/boot.py (не переписываем). Метки: [О] расчёт, [Д] допущение.
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'compression_boot'))
import boot as B
out = {}
def pr(k, v): out[k] = v; print(k, v)

# 0. Воспроизведение базы
r = B.run(15, 'boost', 'base', 1.375)
pr('0_base_Q15_fee1375', dict(capex=r['capex'], npv15=r['npv15'], npv20=r['npv20']))
pr('0_thr_Q15_base', B.solve(lambda x: B.run(15, 'boost', 'base', x)['npv20'], 0.01, 30))
pr('0_thr_Q5_base', B.solve(lambda x: B.run(5, 'boost', 'base', x)['npv20'], 0.01, 30))

# 1. «Ориентир Seplat ×0.6»: установленный капитал 0.6 × пакета ниже цены самого оборудования -> физически невозможно.
#    Что значит $10.8 млн / 50 MMscf/д при разной степени сжатия (вход/выход бар а) [О]:
for ps, pd, st in ((2, 50, 3), (5, 30, 2), (10, 30, 1), (15, 45, 1)):
    hpq = B.hp_per_mmscfd(ps, pd, st); hp_run = hpq*50
    pr(f'1_seplat_{ps}->{pd}bar', dict(hp_per_MMscfd=round(hpq), hp_work=round(hp_run),
       usd_per_hp_work=round(10.8e6/hp_run), x_to_pkg1250=round(10.8e6/hp_run/1250, 2)))

# 2. Оператор «сам» по тем же поправкам аудита (капитал ×2.4-2.8, опекс не ниже скрина ×1.7-2.4).
#    В boot.py оператору дано cm=1.3 и mult=1.0, т.е. аудит снят только с одной стороны сравнения.
def lcoc_oper(Q, conf, sc, r, ov):
    p = dict(B.SC[sc]); o = B.run(Q, conf, sc, operator=True, rate_tax=0.0, ov=ov)
    G = [Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*(p['ramp'] if t == 1 else 1.0) for t in range(1, 11)]
    pvv = sum(g*365/1000/(1+r)**(t+p['build']-1) for t, g in enumerate(G, 1))
    return round(-B.npv(r, o['cf'])/pvv, 3), o['capex']['total'], o['fixed']
for Q in (5, 15):
    for sc in ('opt', 'base'):
        cm = B.SC[sc]['cm']; mult = B.SC[sc]['mult']
        pr(f'2_oper_DIY_Q{Q}_{sc}_as_model(cm1.3)', lcoc_oper(Q, 'boost', sc, .15, {}))
        pr(f'2_oper_DIY_Q{Q}_{sc}_audit(cm{cm},floor x{mult})', lcoc_oper(Q, 'boost', sc, .15, dict(cm=cm, mult=mult)))
        pr(f'2_oper_DIY_Q{Q}_{sc}_audit_r12%', lcoc_oper(Q, 'boost', sc, .12, dict(cm=cm, mult=mult)))

# 3. Объём: станция 3×50% может пропустить до 1.5Q при всплеске дебита. Эмпирика VIIRS (lf_check.py, Q 3-30, годы 1-6):
#    кап 1.0Q: .787 .686 .625 .600 .568 .545; кап 1.5Q: .897 .783 .727 .697 .657 .627 -> ×1.15 к объёму [О].
#    Готовность N+1 по самой станции ~0.97-0.98, но простои приёмного завода/промысла остаются -> 0.95 [Д].
vol = dict(lf=0.78*1.15, up=0.95)
pr('3_thr_Q15_base_vol_fix', B.solve(lambda x: B.run(15, 'boost', 'base', x, ov=vol)['npv20'], 0.01, 30))
pr('3_thr_Q5_base_vol_fix', B.solve(lambda x: B.run(5, 'boost', 'base', x, ov=vol)['npv20'], 0.01, 30))
# + без экспата (местный партнёр) — всё ещё аудит по капиталу ×2.6 и опекс ≥ скрин×2.0
fav = dict(vol, staff=0.25, ovh=0.10)
pr('3_thr_Q15_base_vol_noexpat', B.solve(lambda x: B.run(15, 'boost', 'base', x, ov=fav)['npv20'], 0.01, 30))
pr('3_thr_Q15_base_vol_noexpat_cm2.4', B.solve(lambda x: B.run(15, 'boost', 'base', x, ov=dict(fav, cm=2.4))['npv20'], 0.01, 30))
r = B.run(15, 'boost', 'base', 1.75, ov=fav)
pr('3_Q15_base_vol_noexpat_fee1.75', dict(npv15=r['npv15'], npv20=r['npv20'], e1=r['y1']['ebitda'], e2=r['y2']['ebitda']))

# 4. Аренда + обвязка: бенчмарк USAC в boot.py — только компрессор; оператору нужны сепарация/осушка/учёт/отвод/монтаж.
def bop_usd_per_mscf(Q, sc, cm, r=.15):
    p = B.SC[sc]; bop = 0.25*Q**0.6 + 0.30 + 0.03*Q; cap = bop*cm + 1.0*p['pipe_km'] + 0.15
    crf = r/(1-(1+r)**-10); Gav = sum(Q*p['lf']*(1-p['decl'])**(t-1)*p['up'] for t in range(1, 11))/10*365/1000
    return round((cap*crf + cap*0.05)/Gav, 3)   # + ТОиР/страховка 5% [Д]
for Q in (5, 15):
    for cm in (1.3, 2.6):
        pr(f'4_BoP_Q{Q}_base_cm{cm}', bop_usd_per_mscf(Q, 'base', cm))
pr('4_rental_compressor_only_base', B.contract_compression_benchmark(15, 'boost', 'base', 1.5) , )
pr('4_rental_compressor_only_base_x2.5', B.contract_compression_benchmark(15, 'boost', 'base', 2.5))

# 5. Структура цены: оплата за мощность $/л.с.-мес (как у рынка аренды), без объёмного риска.
def run_cap(Q, sc, rate_hp_mo, hp_basis='work', ov=None):
    p = dict(B.SC[sc]); p.update(ov or {}); cx = B.capex(Q, 'boost', p); C = cx['total']
    hp = cx['hp']/1.5 if hp_basis == 'work' else cx['hp']
    ss = B.staff_scale(Q)
    fixed = max((p['staff']+p['ovh'])*ss + p['sec'] + C*(p['tor']+p['ins']),
                (0.20*ss + 0.10 + cx['screen']*0.025)*p['mult'])
    cf = [-C/p['build']]*p['build']; prev = 0
    for t in range(1, 11):
        G = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*(p['ramp'] if t == 1 else 1.0)
        rev = hp*rate_hp_mo*12/1e6*p['coll']*(1-B.NCDMB)*(p['ramp'] if t == 1 else 1.0)
        e = rev - fixed*1.03**(t-1) - G*365*p['var']/1000
        tx = max(0, (e - C/10)*B.TAX); wc = p['dso']/365*(rev-prev); prev = rev
        cf.append(e - tx - wc + (p['dso']/365*rev if t == 10 else 0))
    return B.npv(.20, cf)
for Q in (5, 15):
    for sc in ('opt', 'base'):
        pr(f'5_capacity_rate_NPV20=0_Q{Q}_{sc}_$per_work_hp_mo', B.solve(lambda x: run_cap(Q, sc, x), 1, 2000))
        pr(f'5_capacity_rate_NPV20=0_Q{Q}_{sc}_noexpat_vol_$per_work_hp_mo',
           B.solve(lambda x: run_cap(Q, sc, x, ov=dict(staff=0.25, ovh=0.10, up=0.95)), 1, 2000))
pr('5_USAC_2025_$per_hp_mo', 21.38)

# 6. Готовность оператора платить при боосте: газ после дожима продаётся оператором -> WTP = штраф/(1-t) + нетбэк газа N.
#    N НЕ НАЙДЕН; [Д] 0-1.0 $/Mscf (DBP ~$2.1-2.4/MMBtu минус переработка и транспорт).
for t in (0.34, 0.51):
    for P in (1.6, 1.9):
        pr(f'6_WTP_t{t}_P{P}', [round(P/(1-t)+N, 2) for N in (0.0, 0.5, 1.0)])
json.dump(out, open(os.path.join(HERE, 'sk_boot_out.json'), 'w'), indent=1, ensure_ascii=False)
