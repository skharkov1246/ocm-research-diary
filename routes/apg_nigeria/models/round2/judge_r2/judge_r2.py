# Раунд 2, судья: числа, которых не хватало в решении совета (критика, разделы (б) и (в)).
# Всё считается на исправленной модели сделки ../deal/deal_r2.py (4.40 МВт на MMscf/д, загрузка 0.78, падение 8%/г,
# курсовой лаг, DoP возвращает оплату газа + эскроу, налог 34%). Запуск из любого каталога:
#   python3 routes/apg_nigeria/models/round2/judge_r2/judge_r2.py  -> печать + judge_r2_out.json рядом
# Метки: [О] -- расчёт этого скрипта; [Д] -- допущение (диапазон в комментарии); [В]/[П] -- источник.
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'deal'))
from deal_r2 import NEW, STRUCTS, TL_NEW, run, be, site, station, FX, npv

S = STRUCTS; TL = TL_NEW; OUT = {}
OPF = 0.15/1.05                       # плата оператору $0.15/Mscf в $/MMBtu [Д, НЕ НАЙДЕНО]
def g_allin(a): return a-OPF          # цена газа «всё включено» -> параметр gas модели
def kw_e(e): return dict(edti=e, depyrs=5 if e else 10)
B = {'1y': (dict(TL, delay=1), {}), '2y': (dict(TL, delay=2), dict(capex_split=[0.5, 0.5])),
     '2y_slip': (dict(TL, delay=2), {})}
def R(T, e, b, gas=0.74, sn='25km', s='dopgas', **kw):
    tl, ex = B[b]
    return run(dict(site(sn), p_dir=T, gas=g_allin(gas)), S[s], tl, **ex, **kw_e(e), **kw)

# ------------------------------------------------------------------ 1. Ожидаемая стоимость программы ворот
# Вероятности [Д] -- те же, что в ../dev_sale/dev_sale.py (данных нет; у 9 факелов 0 подписанных покупателей):
#   p0 -- LOI якоря по цене >= T_LOI; p1 -- опцион/СП с держателем NGFCP или GCA с NEPL + согласие NUPRC;
#   p2 -- 6 мес. замера >= порога и подписанный PPA; pfc -- финансовое закрытие (кредитор/DFI).
SC = {'pess': dict(p0=0.10, p1=0.40, p2=0.35, pfc=0.50, k=1.0),
      'base': dict(p0=0.20, p1=0.55, p2=0.50, pfc=0.65, k=0.5),
      'opt':  dict(p0=0.30, p1=0.70, p2=0.65, pfc=0.80, k=0.0)}
# Внешние статьи ворот [О: сумма статей решения судьи, Д по статьям]; команда, поездки, охрана полевых работ
# $15-40 тыс./мес [Д, как в dev_sale.py] по длительности ворот нового графика: 0 -- 2 мес., 1 -- 4 мес., 2 -- 8 мес.
# ворота 0: обход $20-50k + юрист $15-40k + данные NGFCP $2k обзор + $1k за площадку только в программе ($2-11k)
#           + котировки СПГ/CNG $0-10k + MAIN $0.5k = $38-112k [О: сумма статей]
EXT = dict(g0=(0.038, 0.112), g1=(0.110, 0.270), g2=(0.272, 0.442))
TEAM = dict(g0=(0.030, 0.080), g1=(0.060, 0.160), g2=(0.120, 0.320))
TT = dict(g0=0.10, g1=0.35, g2=0.85)   # середина трат, лет от 09.10.2026
T_FID = 1.2                            # FID 12.2027 [график раунда 2]
def cost(g, k): return (EXT[g][0]+TEAM[g][0])*(1-k)+(EXT[g][1]+TEAM[g][1])*k
def ev(npv_fid, sc, r):
    pr = dict(g0=1.0, g1=sc['p0'], g2=sc['p0']*sc['p1'])
    c = sum(pr[g]*cost(g, sc['k'])/(1+r)**TT[g] for g in ('g0', 'g1', 'g2'))
    pf = sc['p0']*sc['p1']*sc['p2']*sc['pfc']
    return pf*npv_fid/(1+r)**T_FID - c, pf, c
prog = {}
for label, (T, e, b) in {'LOI ₦495 без EDTI, 2 г': (495, False, '2y'),
                          'нижняя ступень ₦440, EDTI, 2 г': (440, True, '2y'),
                          '₦405 EDTI 1 г (лучший график)': (405, True, '1y'),
                          '₦550 без EDTI, 2 г': (550, False, '2y'),
                          '₦600 без EDTI, 2 г': (600, False, '2y')}.items():
    for gas in (0.74, 1.5):
        x = R(T, e, b, gas)
        for scn, sc in SC.items():
            e15, pf, c15 = ev(x['npv15'], sc, 0.15); e20, _, c20 = ev(x['npv20'], sc, 0.20)
            prog[f'{label}|газ {gas}|{scn}'] = dict(npv15_fid=x['npv15'], npv20_fid=x['npv20'], p_fid=round(pf, 4),
                                                   cost_pv15=round(c15, 3), ev15=round(e15, 3), ev20=round(e20, 3))
OUT['program_ev'] = prog
# безубыточная p0 (вероятность LOI) при прочих базовых вероятностях
def p0_star(npv_fid, r):
    lo, hi = 0.0, 1.0
    f = lambda p: ev(npv_fid, dict(SC['base'], p0=p), r)[0]
    if f(hi) < 0: return None
    for _ in range(60):
        m = (lo+hi)/2
        if f(m) < 0: lo = m
        else: hi = m
    return round(hi, 3)
x495 = R(495, False, '2y'); x440 = R(440, True, '2y')
OUT['p0_breakeven'] = {'₦495 без EDTI 2 г, газ 0.74': dict(ev15=p0_star(x495['npv15'], .15), ev20=p0_star(x495['npv20'], .20)),
                       '₦440 EDTI 2 г, газ 0.74': dict(ev15=p0_star(x440['npv15'], .15), ev20=p0_star(x440['npv20'], .20))}
# тариф, при котором программа даёт EV20 = 0 при базовых вероятностях (без EDTI, 2 г, газ 0.74)
lo, hi = 400.0, 1200.0
for _ in range(40):
    m = (lo+hi)/2
    if ev(R(m, False, '2y')['npv20'], SC['base'], .20)[0] < 0: lo = m
    else: hi = m
OUT['T_program_ev20_zero_base'] = round(hi)
print('=== 1. Программа ворот: P(FID) и ожидаемая стоимость, $ млн [О; вероятности Д]')
for k, v in prog.items():
    if 'газ 0.74' in k: print('  ', k, v)
print('   безубыточная P(LOI) при прочих базовых:', OUT['p0_breakeven'])
print('   тариф, при котором EV20 программы = 0 (база, без EDTI, 2 г):', OUT['T_program_ev20_zero_base'])

# ------------------------------------------------------------------ 2. Бюджет до FID с командой
team_total = (sum(TEAM[g][0] for g in TEAM), sum(TEAM[g][1] for g in TEAM))
ext_total = (sum(EXT[g][0] for g in EXT), sum(EXT[g][1] for g in EXT))
OUT['budget_to_fid'] = dict(external=ext_total, team=team_total,
                            total=(round(ext_total[0]+team_total[0], 3), round(ext_total[1]+team_total[1], 3)),
                            gate0_with_team=(round(EXT['g0'][0]+TEAM['g0'][0], 3), round(EXT['g0'][1]+TEAM['g0'][1], 3)))
print('\n=== 2. До FID, $ млн:', OUT['budget_to_fid'])

# ------------------------------------------------------------------ 3. Цена норм регуляторного пакета на модели раунда 2
norms = {}
# Норма 1: тариф в найре без индексации против долларового, при T_LOI (NPV15, $ млн); матрица §11 считает это для 1 г;
# здесь -- для ячейки ворот 0 (₦495, без EDTI, плановая стройка 2 г).
base495 = R(495, False, '2y')
for dep in (0.08, 0.12, 0.16):
    tl, ex = B['2y']
    tln = dict(tl, idx=0.0, dep=dep, fx_lag_q=0.0, fx_pay=False)
    y = run(dict(site('25km'), p_dir=495, gas=g_allin(0.74)), S['dopgas'], tln, **ex, **kw_e(False))
    norms[f'N1|₦495 2г без EDTI|найра, девальвация {int(dep*100)}%'] = dict(npv15=y['npv15'], npv20=y['npv20'],
        d_npv15=round(base495['npv15']-y['npv15'], 2))
# Норма 2: DoP в форме компенсации маржи против DoP «возврат оплаты газа» (ячейка ворот 0)
dm = R(495, False, '2y', s='dopmargin')
norms['N2|₦495 2г без EDTI|DoP-маржа против DoP-газ'] = dict(d_npv15=round(dm['npv15']-base495['npv15'], 2),
    d_npv20=round(dm['npv20']-base495['npv20'], 2),
    d_be20=round(be(dict(site('25km'), gas=g_allin(0.74)), S['dopgas'], B['2y'][0], 'npv20', **B['2y'][1], **kw_e(False))
                 - be(dict(site('25km'), gas=g_allin(0.74)), S['dopmargin'], B['2y'][0], 'npv20', **B['2y'][1], **kw_e(False))))
# Норма 3: EDTI (та же ячейка)
ed = R(495, True, '2y')
norms['N3|₦495 2г|EDTI против без EDTI'] = dict(d_npv15=round(ed['npv15']-base495['npv15'], 2), d_npv20=round(ed['npv20']-base495['npv20'], 2))
# Норма 4: часть штрафа платит утилизатору: $1.0-1.75/Mscf x 0.6 [Д: коэффициент реализации раунда 1] на сожжённый газ с 1-го года
for fee in (1.0, 1.75):
    for (lab, e, b) in (('ворота 0: без EDTI 2 г', False, '2y'), ('EDTI 1 г', True, '1y')):
        tl, ex = B[b]
        p = dict(site('25km'), gas=g_allin(0.74))
        b0 = be(p, S['dopgas'], tl, 'npv20', **ex, **kw_e(e))
        b1 = be(p, S['dopgas'], tl, 'npv20', earnout=(-0.6*fee, 1), **ex, **kw_e(e))
        norms[f'N4|{lab}|плата ${fee}/Mscf x 0.6'] = dict(be20_before=round(b0), be20_after=round(b1), d_be20=round(b1-b0))
# Норма 5: порог eligible customer 1-2 МВт·ч/ч. Площадка 1 MMscf/д: 1 x 0.78 x 4.40 x 0.96 = 3.29 МВт нетто.
q1 = {}
for lab, km, dev in (('трасса 10 км, общая часть $2.0 млн', 10.0, 0.0), ('трасса 25 км, общая часть $2.8 млн', 25.0, 0.8)):
    for e, b in ((True, '1y'), (False, '2y')):
        tl, ex = B[b]
        p = dict(site('25km'), Q=1.0, anchor=3.29, line_km=km, loss=0.03*km/25, dev=dev, gas=g_allin(0.74))
        q1[f'{lab}|{"EDTI 1г" if e else "без EDTI 2г"}'] = dict(capex=station(p)['C'].__round__(2),
            be20=round(be(p, S['dopgas'], tl, 'npv20', hi=4000.0, **ex, **kw_e(e))))
        pl = dict(p, staff_mult=1.0)   # «бережливо»: персонал без экспата, охрана и община как есть [Д]
        q1[f'{lab}|{"EDTI 1г" if e else "без EDTI 2г"}|персонал x1.0'] = dict(
            be20=round(be(pl, S['dopgas'], tl, 'npv20', hi=4000.0, **ex, **kw_e(e))))
norms['N5|1 MMscf/д, 3.29 МВт'] = q1
OUT['norms'] = norms
print('\n=== 3. Цена норм на модели раунда 2 [О]')
for k, v in norms.items(): print('  ', k, v)

# ------------------------------------------------------------------ 4. Цена раскрытия порогов (что теряем в переговорах)
disc = {}
x74 = R(495, False, '2y', 0.74); x15 = R(495, False, '2y', 1.5)
disc['газ по нашему стопу $1.5 вместо цели $0.74 (₦495, 2 г, без EDTI)'] = dict(d_npv15=round(x15['npv15']-x74['npv15'], 2),
                                                                                d_npv20=round(x15['npv20']-x74['npv20'], 2))
x505 = R(505, False, '2y', 0.74)
disc['каждые ₦10 тарифа выше T_LOI'] = dict(d_npv15=round(x505['npv15']-x74['npv15'], 2), d_npv20=round(x505['npv20']-x74['npv20'], 2))
# потолок earn-out при ₦495 и газе $0.74 (путь А) и потеря, если держатель забирает его целиком
lo, hi = 0.0, 20.0
for _ in range(50):
    m = (lo+hi)/2
    if R(495, False, '2y', 0.74, earnout=(m, 2))['npv20'] >= 0: lo = m
    else: hi = m
eo_cap = round(lo, 2)
xe03 = R(495, False, '2y', 0.74, earnout=(0.30, 2)); xec = R(495, False, '2y', 0.74, earnout=(eo_cap, 2))
disc['earn-out держателю: потолок вместо $0.30/Mscf'] = dict(cap_usd_mscf=eo_cap, d_npv15=round(xec['npv15']-xe03['npv15'], 2),
                                                             d_npv20=round(xec['npv20']-xe03['npv20'], 2))
# EPC видит потолок капитала: капитал при NPV20=0 (газ 0.74) вместо базы x2.6
lo, hi = 2.0, 5.0
for _ in range(40):
    m = (lo+hi)/2
    tl, ex = B['2y']
    if run(dict(site('25km'), p_dir=495, gas=g_allin(0.74), cm=m), S['dopgas'], tl, **ex, **kw_e(False))['npv20'] >= 0: lo = m
    else: hi = m
xc = run(dict(site('25km'), p_dir=495, gas=g_allin(0.74), cm=lo), S['dopgas'], B['2y'][0], **B['2y'][1], **kw_e(False))
disc['EPC-цена по нашему потолку капитала'] = dict(cm=round(lo, 2), capex=xc['capex'], d_npv15=round(xc['npv15']-x74['npv15'], 2),
                                                    d_npv20=round(xc['npv20']-x74['npv20'], 2))
disc['запас NPV20 при ₦495 и газе $0.74 (весь избыток над 20%)'] = dict(npv15=x74['npv15'], npv20=x74['npv20'])
OUT['disclosure'] = disc
print('\n=== 4. Цена раскрытия (₦495, без EDTI, 2 г) [О]')
for k, v in disc.items(): print('  ', k, v)

# ------------------------------------------------------------------ 5. Охрана: блокада/захват на 3 мес. -- допустимая частота
loss20 = []; loss15 = []
for t in range(1, 11):
    y = R(495, False, '2y', 0.74, ev={t: 0.75})
    loss20.append(x74['npv20']-y['npv20']); loss15.append(x74['npv15']-y['npv15'])
h20 = x74['npv20']/sum(loss20)        # годовая вероятность 3-мес. блокады, при которой ожидаемое NPV20 = 0 (малые h)
OUT['blockade'] = dict(loss_npv20_by_year=[round(v, 2) for v in loss20], loss_npv15_by_year=[round(v, 2) for v in loss15],
                       h_max_per_year_npv20=round(h20, 3),
                       note='газ $0.74; при газе на стопе запас NPV20 ~0, и допустимая частота ~0')
print('\n=== 5. Блокада 3 мес. в год t (₦495, без EDTI, 2 г): потеря NPV20', OUT['blockade']['loss_npv20_by_year'],
      '| допустимая частота', OUT['blockade']['h_max_per_year_npv20'], 'в год')

# ------------------------------------------------------------------ 6. Счёт покупателя, аккредитив, проценты на стройке
x = R(495, False, '2y', 0.74)
bill = {gwh: round(gwh*1e6*495/FX/1e6, 2) for gwh in (56.8, 61.9)}
OUT['buyer_bill_495'] = dict(annual_musd=bill, lc_3m_musd={k: round(v/4, 2) for k, v in bill.items()})
C = x['capex']; D = 0.6*C
idc_1y = D*0.14*0.5; idc_2y = (0.5*D/2)*0.14 + (0.5*D+0.5*D/2)*0.14
OUT['idc_rough'] = dict(capex=C, debt=round(D, 2), idc_1y=round(idc_1y, 2), idc_2y_50_50=round(idc_2y, 2),
                        to_finance_1y=round(C+1.3+idc_1y+2.0, 1), to_finance_2y=round(C+1.3+idc_2y+2.0, 1),
                        note='оборотный капитал $1.3 млн и DSRA $2.0 млн -- из раунда 1 [О раунда 1]; грубо')
print('\n=== 6. Счёт покупателя при ₦495, $ млн/г:', OUT['buyer_bill_495'], '| к финансированию:', OUT['idc_rough'])

# ------------------------------------------------------------------ 7. Норма 5: экономический пол якоря без юридического
# МВт нетто, при котором BE20 (газ на стопе $1.5) <= 0.85 x альтернатива покупателя [Д: экономия 15%]; 25 км, DoP-газ
eco = {}
for lab, e, b in (('ворота 0: без EDTI, 2 г', False, '2y'), ('EDTI, 1 г', True, '1y')):
    tl, ex = B[b]
    for pd in (530, 650, 770):
        f = lambda mw: be(dict(site('25km'), anchor=mw, gas=g_allin(1.5)), S['dopgas'], tl, 'npv20', **ex, **kw_e(e)) - 0.85*pd
        if f(15.0) > 0: eco[f'{lab}|дизель {pd}'] = None; continue
        lo, hi = 2.0, 15.0
        if f(lo) <= 0: eco[f'{lab}|дизель {pd}'] = '<2'; continue
        for _ in range(30):
            m = (lo+hi)/2
            if f(m) > 0: lo = m
            else: hi = m
        eco[f'{lab}|дизель {pd}'] = round(hi, 1)
OUT['norm5_econ_floor_mw'] = eco
print('\n=== 7. Экономический пол якоря (МВт нетто, газ на стопе, экономия 15%), 25 км:', eco)

json.dump(OUT, open(os.path.join(HERE, 'judge_r2_out.json'), 'w'), ensure_ascii=False, indent=1)
print('\nзаписано', os.path.join(HERE, 'judge_r2_out.json'))
