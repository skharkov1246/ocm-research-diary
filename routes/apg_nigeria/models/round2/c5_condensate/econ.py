# Раунд 2. Экономика «Стабилизация конденсата C5+ (JT-сепарация) с подачей в нефтяной поток оператора».
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
# Поправки аудита NaCN: капитал x2.4/2.6/2.8 к пакету оборудования; денежный опекс >= скрин x1.7/2.0/2.4;
# персонал с экспатом; ТОиР 3.5-5%; выход по нижней границе; спрос /1.3-2.7 до контракта (здесь — доля
# газа высокого давления и доля сделки до подписания); второй год: стройка, разгон, загрузка 0.78 (VIIRS),
# падение дебита 8%/г, оборотный капитал, девальвация на дебиторке при оплате в найре, налог 34%.
# Выходы по составу берутся из process_out.json (PR-флэш, process.py).
# Запуск: python3 econ.py -> econ_out.json и econ_stdout.txt рядом.
import json, math, os, csv
HERE = os.path.dirname(os.path.abspath(__file__))
TAX = 0.34                      # CIT 30% + development levy 4% (NTA 2025) [В]
YEARS = 10                      # лет эксплуатации [Д 8-12]
PROC = json.load(open(os.path.join(HERE, 'process_out.json')))
NOINH = json.load(open(os.path.join(HERE, 'noinhib_out.json')))

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))
def irr(cf):
    lo, hi = -0.99, 3.0
    if npv(lo, cf)*npv(hi, cf) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if npv(lo, cf)*npv(m, cf) <= 0: hi = m
        else: lo = m
    return m

# ---------- выход, барр. стабильного конденсата на 1 MMscf газа через скид ----------
def jt_yield(comp, key):
    """key: 'P1->P2' — глубокий JT с рекуперацией и ингибитором; 'noinh:P1' — холодный сепаратор ВД при Th+3 C без ингибитора"""
    if key.startswith('noinh:'):
        return NOINH[comp][key.split(':')[1]]['yield_stab_bbl_per_MMscf']
    r = PROC['jt_recup'].get(comp, {}).get(key)
    return 0.0 if not r else r['yield_stab_bbl_per_MMscf']
def meoh_kg(comp, key, hi=False):
    if key.startswith('noinh:'): return 0.0
    r = PROC['jt_recup'].get(comp, {}).get(key)
    if not r or not r['meoh'].get('need'): return 0.0
    return r['meoh']['kg_per_MMscf_hi' if hi else 'kg_per_MMscf_lo']
# классы газа: средний = пробы SPE AG2/AG6, жирный = SPE AG4 / NJTD, тощий = lean_ND (и SPE AG1) [П/В].
# База — нижняя граница: проба с меньшим C5+ и схема без ингибитора (реагентов нет) при газе ВД 30 бар.
# Оптимизм — глубокий JT 50->3 бар с рекуперацией и метанолом; пессимизм — без ингибитора при 15 бар.
GAS = {
 'mid':  dict(opt=('SPE_AG6', '50->3'), base=('SPE_AG2', 'noinh:30'), pess=('SPE_AG2', 'noinh:15')),
 'rich': dict(opt=('NJTD_A', '50->3'),  base=('SPE_AG4', 'noinh:30'), pess=('SPE_AG4', 'noinh:15')),
 'lean': dict(opt=('lean_ND', '50->3'), base=('lean_ND', 'noinh:30'), pess=('lean_ND', 'noinh:15')),
}

SC = {
 'opt':  dict(lf=0.85, decl=0.05, up=0.97, hp=1.0, build=1, ramp=0.8, cm=2.4, dev=0.30, opx_mult=1.7,
              brent=85, diff=0.0, roy=0.05, loss=0.03, tariff=1.5, gas_pay='shrink', gas_usd_mmbtu=0.5,
              meoh_usd_t=600, meoh_hi=False, staff=0.30, sec=0.15, comm=0.05, tor=0.035, ins=0.010, ovh=0.10,
              dso=45, coll=0.98, fxloss=0.0, deal=1.0),
 'base': dict(lf=0.78, decl=0.08, up=0.95, hp=0.6, build=1, ramp=0.7, cm=2.6, dev=0.50, opx_mult=2.0,
              brent=70, diff=3.0, roy=0.17, loss=0.08, tariff=3.0, gas_pay='shrink', gas_usd_mmbtu=1.0,
              meoh_usd_t=750, meoh_hi=False, staff=0.40, sec=0.25, comm=0.08, tor=0.0425, ins=0.015, ovh=0.15,
              dso=75, coll=0.95, fxloss=0.0, deal=1/2.0),
 'pess': dict(lf=0.65, decl=0.12, up=0.92, hp=0.3, build=2, ramp=0.5, cm=2.8, dev=0.80, opx_mult=2.4,
              brent=55, diff=8.0, roy=0.185, loss=0.15, tariff=5.0, gas_pay='throughput', gas_usd_mmbtu=0.25,
              meoh_usd_t=900, meoh_hi=True, staff=0.50, sec=0.40, comm=0.12, tor=0.05, ins=0.020, ovh=0.20,
              dso=120, coll=0.85, fxloss=0.12, deal=1/2.7),
}
# Пояснения к допущениям [Д]:
#  lf — загрузка факела по VIIRS (0.78 по корпусу), decl — падение дебита (медиана VIIRS 6 лет -9.2%/г [О]);
#  hp — доля факельного газа, доступного на входе скида при давлении P1 (газ сепараторов ВД). НЕ НАЙДЕНО по площадкам.
#  brent — долгосрочная цена Brent, $/барр.; spot 10.2026 ~$90-103 [В], прогноз Fitch на 2027 $71 [В];
#  roy — роялти, если конденсат считается нефтью в руках лицензиата: 15% суша + ценовое роялти ~2% при $70
#        [В: EY 2021, Deloitte]; оптимизм — режим NGL 5% [В]; пессимизм — 15% + 3.5%.
#  loss — потери в магистрали (Trans Forcados 2015: 10-21% закачки, Mart Resources [В]; 2025 ниже [В NUPRC]);
#  tariff — плата оператору за приём/транспорт, $/барр. НЕ НАЙДЕНО (договоры crude handling не публичны).
#  gas_pay — 'shrink': платим только за энергию, ушедшую в жидкость; 'throughput': $/Mscf на весь газ через скид
#        (пол NGFCP $0.25/Mscf [В корпус]).
#  deal — доля вероятности/объёма до подписанного договора с оператором (спрос /1.3-2.7).
MMBTU_PER_BBL_C5 = 4.6          # HHV стабильного C5+ ~4.6-5.0 MMBtu/барр. [О: 4380 Btu/scf C5+ пара, ~5.2 бар газа/галлон]

def package(Q, lp_compression=False, p_dis=30.0):
    s = (Q/5.0)**0.6
    pk = dict(jt_lts_hx=0.45*s,          # входной сепаратор, газ/газ теплообменник, JT-клапан, НТС, впрыск ингибитора [Д 0.3-0.7]
              stabilizer=0.30*s,         # подогреваемый стабилизатор/колонна с ребойлером на топливном газе [Д 0.2-0.5]
              tank_pumps_line=0.25,      # буферный резервуар ~3 сут, насосы, врезка в нефтяной коллектор ≤1 км [Д 0.15-0.4]
              metering=0.30,             # учётный узел (Кориолис, пробоотборник) под требования NUPRC [Д 0.2-0.5]
              ic_tieins=0.15*s)          # КИПиА, врезки в факельный коллектор [Д]
    if lp_compression:
        R = (p_dis/2.0)**(1/3); hpq = 22*R*3*1.10   # GPSA-эмпирика, как в g2p.py [В/О]
        pk['compression'] = hpq*Q*1.5*1250/1e6      # N+1, $1250/л.с. пакет [Д 1000-1500]
    return pk

def run(Q, gas='mid', sc='base', ov=None, years=YEARS, contract=True, lp_compression=False, operator=False):
    p = dict(SC[sc]); p.update(ov or {})
    comp, key = GAS[gas][sc]
    comp = p.get('comp', comp); key = p.get('jt_key', key)
    y = p.get('yield', jt_yield(comp, key))                         # барр./MMscf
    mk = p.get('meoh_kg', meoh_kg(comp, key, p['meoh_hi']))         # кг MeOH/MMscf
    pk = package(Q, lp_compression); pks = sum(pk.values())
    cm = p['cm'] if not operator else 1.3
    dev = p['dev'] if not operator else 0.15
    inst = pks*cm + dev
    idc = inst*0.14*p['build']/2
    C = inst + idc
    staff_sc = 1.0 if Q <= 1.5 else (1.25 if Q <= 7 else 1.5)
    if operator:   # оператор делает сам: персонал, охрана и площадка уже есть [Д]
        fixed_bu = 0.08*staff_sc + 0.05 + C*(0.035+0.01)
        fixed_scr = fixed_bu
    else:
        fixed_bu = (p['staff']+p['ovh'])*staff_sc + p['sec'] + p['comm'] + C*(p['tor']+p['ins'])
        fixed_scr = 0.12*staff_sc + 0.05 + pks*(0.02+0.005)          # промоутерский скрин [Д]
    fixed = max(fixed_bu, fixed_scr*(p['opx_mult'] if not operator else 1.0))
    if lp_compression and not operator: fixed += 0.10*staff_sc       # механик по компрессорам [Д]
    b = p['build']; cf = [-C/b]*b; rows = []; prev = 0.0
    dem = 1.0 if contract else p['deal']
    net_px = (p['brent'] - p['diff'])*(1 - p['roy'])*(1 - p['loss']) - p['tariff']   # $/барр. у нас
    for t in range(1, years+1):
        G = Q*p['lf']*(1-p['decl'])**(t-1)*p['up']*p['hp']*(p['ramp'] if t == 1 else 1.0)*dem   # MMscf/д через скид
        bbl = G*y*365
        rev = bbl*net_px/1e6*p['coll']*(1 - p['fxloss']*p['dso']/365)
        if p['gas_pay'] == 'shrink': gas_c = bbl*MMBTU_PER_BBL_C5*p['gas_usd_mmbtu']/1e6
        else: gas_c = G*1000*365*p['gas_usd_mmbtu']/1e6
        meoh = G*365*mk/1000*p['meoh_usd_t']/1e6
        comp_fuel_om = (G*1000*365*0.055*0.0/1e6) if lp_compression else 0.0           # топливо из факела бесплатно [Д]
        opex = fixed*(1.03)**(t-1) + gas_c + meoh + comp_fuel_om
        e = rev - opex
        tx = max(0.0, (e - C/10)*TAX)
        wc = p['dso']/365*(rev-prev); prev = rev
        c = e - tx - wc + (p['dso']/365*rev if t == years else 0)
        cf.append(c)
        rows.append(dict(t=t, G=round(G, 2), bbl_d=round(G*y, 1), rev=round(rev, 3), opex=round(opex, 3), ebitda=round(e, 3)))
    i = irr(cf); acc = 0; pb = None
    for k, c in enumerate(cf):
        if pb is None and acc < 0 and acc + c >= 0 and k > 0: pb = round(k-1 + (-acc)/c - (b-1), 1)
        acc += c
    return dict(Q=Q, gas=gas, sc=sc, comp=comp, jt=key, yield_bbl_per_MMscf=y, meoh_kg_per_MMscf=mk,
                net_usd_bbl=round(net_px, 2), package=round(pks, 2), capex=round(C, 2), fixed=round(fixed, 3),
                y1=rows[0], y2=rows[1], y5=rows[4],
                npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2),
                irr=None if i is None else round(100*i, 1), payback=pb)

def solve(fn, lo, hi, it=70):
    flo, fhi = fn(lo), fn(hi)
    if (flo > 0) == (fhi > 0): return None
    for _ in range(it):
        m = (lo+hi)/2; fm = fn(m)
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
    return round((lo+hi)/2, 2)

# Именованные факелы (VIIRS 2025 / минимум 2023-25, MMscf/д) [О по wb/nigeria_flares.csv, World Bank П]
def viirs():
    H = {}; M = {}
    for r in csv.DictReader(open(os.path.join(HERE, '..', '..', 'wb', 'nigeria_flares.csv'))):
        try: q = float(r['MMSCFD']); yy = int(r['Year'])
        except ValueError: continue
        H.setdefault(r['id'], {}); H[r['id']][yy] = H[r['id']].get(yy, 0) + q
        M[r['id']] = (r['Field name'], r['Operator'], r['Location'])
    return H, M
NAMED = ['Umuseti', 'Utorogu', 'Ughelli East', 'Sapele', 'Oben', 'Kwale', 'Opuama', 'Oredo', 'Obiafu-Obrikom']

if __name__ == '__main__':
    out = dict(generic={}, named={}, thresholds={}, addon={}, operator={}, lp={}, sens={})
    for Q in (1, 5, 15):
        for gas in ('lean', 'mid', 'rich'):
            for sc in ('opt', 'base', 'pess'):
                out['generic'][f'Q{Q}_{gas}_{sc}'] = run(Q, gas, sc)
            out['generic'][f'Q{Q}_{gas}_base_precontract'] = run(Q, gas, 'base', contract=False)
            out['operator'][f'Q{Q}_{gas}_base'] = run(Q, gas, 'base', operator=True)
            out['lp'][f'Q{Q}_{gas}_base'] = run(Q, gas, 'base', lp_compression=True, ov=dict(hp=1.0))
    H, M = viirs()
    for name in NAMED:
        ids = [k for k, v in M.items() if v[0] == name and v[2] == 'ONSHORE']
        q25 = sum(H[i].get(2025, 0) for i in ids); qmin = min(sum(H[i].get(yy, 0) for i in ids) for yy in (2023, 2024, 2025))
        op = M[ids[0]][1] if ids else '?'
        for qn, q in (('q2025', q25), ('qmin', qmin)):
            for gas in ('mid', 'rich'):
                for sc in ('opt', 'base', 'pess'):
                    r = run(max(q, 0.05), gas, sc)
                    out['named'][f'{name} ({op})|{qn}={q:.1f}|{gas}|{sc}'] = {k: r[k] for k in ('capex', 'y2', 'npv15', 'npv20', 'irr', 'payback')}
    # пороги NPV20 = 0 (база, договор подписан)
    for Q in (1, 5, 15):
        for gas in ('mid', 'rich'):
            for sc in ('base', 'opt'):
                out['thresholds'][f'Q{Q}_{gas}_{sc}'] = dict(
                    brent=solve(lambda x: run(Q, gas, sc, ov=dict(brent=x))['npv20'], 20, 2000),
                    yield_bbl_per_MMscf=solve(lambda x: run(Q, gas, sc, ov=dict(yield_=0, **{'yield': x}))['npv20'], 0.01, 200),
                    capex_mult=solve(lambda x: run(Q, gas, sc, ov=dict(cm=x))['npv20'], 0.05, 5.0),
                    Q_mmscfd=solve(lambda x: run(x, gas, sc)['npv20'], 0.2, 400))
    # надстройка к «энергии якорю»: газ станции 11/(4.40*0.96)=2.60 MMscf/д идёт через узел точки росы,
    # который нужен двигателям и так (капитал JT-части уже в $2.8 млн подготовки газа). Инкремент: стабилизатор,
    # резервуар, учёт, врезка. Персонал и охрана общие -> инкремент 0.10-0.15 $M/г [Д].
    MMSCFD_STATION = 11/(4.40*0.96)
    for gas in ('mid', 'rich'):
        for sc in ('opt', 'base', 'pess'):
            p = dict(SC[sc]); comp, key = GAS[gas][sc]; y = jt_yield(comp, key)
            if key.startswith('noinh'):  # газ станции подаётся к двигателям под давлением: берём глубину схемы двигателя = без ингибитора
                pass
            pk = 0.30*(MMSCFD_STATION/5)**0.6 + 0.25 + 0.30           # стабилизатор + резервуар/линия + учёт [Д]
            C = pk*p['cm'] + 0.15
            net_px = (p['brent']-p['diff'])*(1-p['roy'])*(1-p['loss']) - p['tariff']
            cf = [-C]
            for t in range(1, YEARS+1):
                G = MMSCFD_STATION*0.92*(1-p['decl'])**(t-1)*(0.8 if t == 1 else 1)   # готовность станции 0.92 [deal]
                rev = G*y*365*net_px/1e6*p['coll']
                e = rev - (0.12 + C*(p['tor']+p['ins']))
                cf.append(e - max(0, (e-C/10)*TAX))
            out['addon'][f'power_anchor_{gas}_{sc}'] = dict(station_mmscfd=round(MMSCFD_STATION, 2), yield_=y, bbl_d=round(MMSCFD_STATION*0.92*y, 1),
                capex=round(C, 2), ebitda_y2=round(cf[2] + 0, 3), npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2))
    # чувствительности, Q=5 средний газ база
    base = dict(Q=5, gas='mid', sc='base')
    for lab, ov in (('base', {}), ('brent_90', dict(brent=90)), ('brent_100', dict(brent=100)), ('roy_5', dict(roy=0.05)),
                    ('loss_0_tariff_0', dict(loss=0.0, tariff=0.0)), ('hp_1.0', dict(hp=1.0)), ('cm_1.0', dict(cm=1.0)),
                    ('gas_throughput_0.25', dict(gas_pay='throughput', gas_usd_mmbtu=0.25)), ('no_meoh', dict(meoh_kg=0.0)),
                    ('staff_sec_half', dict(staff=0.2, sec=0.12)), ('deep_jt_30->3_meoh', dict(jt_key='30->3')), ('deep_jt_50->3_meoh', dict(jt_key='50->3')), ('noinh_15', dict(jt_key='noinh:15')), ('noinh_50', dict(jt_key='noinh:50')), ('AG6_noinh30', dict(comp='SPE_AG6')), ('AG6_deep30_meoh', dict(comp='SPE_AG6', jt_key='30->3'))):
        r = run(5, 'mid', 'base', ov=ov)
        out['sens'][lab] = {k: r[k] for k in ('yield_bbl_per_MMscf', 'net_usd_bbl', 'capex', 'fixed', 'y2', 'npv15', 'npv20', 'irr')}
    json.dump(out, open(os.path.join(HERE, 'econ_out.json'), 'w'), ensure_ascii=False, indent=1)
    lines = []
    for sec in ('generic', 'operator', 'lp', 'addon', 'thresholds', 'sens'):
        lines.append(f'== {sec}')
        for k, v in out[sec].items():
            if isinstance(v, dict) and 'y2' in v and sec != 'sens':
                lines.append(f"{k}: y={v['yield_bbl_per_MMscf']} net${v['net_usd_bbl']} pkg{v['package']} capex{v['capex']} fixed{v['fixed']} "
                             f"y2:G{v['y2']['G']} bbl/d{v['y2']['bbl_d']} rev{v['y2']['rev']} opex{v['y2']['opex']} E{v['y2']['ebitda']} "
                             f"NPV15 {v['npv15']} NPV20 {v['npv20']} IRR {v['irr']} pb {v['payback']}")
            else: lines.append(f'{k}: {v}')
    lines.append('== named')
    for k, v in out['named'].items(): lines.append(f'{k}: {v}')
    open(os.path.join(HERE, 'econ_stdout.txt'), 'w').write('\n'.join(lines))
    print('\n'.join(lines))
