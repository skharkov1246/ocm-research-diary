# Раунд 2. Конденсат C5+ с JT-сепарации на факельном газе Нигерии: выход по составу газа.
# Метки: [П] первичный, [В] вторичный, [О] расчёт, [Д] допущение с диапазоном.
# Уравнение состояния Пенга-Робинсона (пакет thermo 0.6.1, kij ChemSep PR) [О].
# Запуск: <venv с thermo>/bin/python process.py  -> process_out.json рядом.
import json, math, os
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
HERE = os.path.dirname(os.path.abspath(__file__))
KEYS = ['C1','C2','C3','iC4','nC4','iC5','nC5','C6','C7','C8','C9','C10','N2','CO2']
IDS  = ['methane','ethane','propane','isobutane','butane','isopentane','pentane','hexane','heptane',
        'octane','nonane','decane','nitrogen','carbon dioxide']
# Плотность идеальной жидкости, отн. к воде при 60F (GPA 2145 / таблица SPE-184345) [П]
SG = dict(C1=0.300,C2=0.356,C3=0.507,iC4=0.563,nC4=0.584,iC5=0.625,nC5=0.631,C6=0.690,C7=0.727,C8=0.749,C9=0.768,C10=0.782,N2=0.809,CO2=0.818)
KMOL_PER_MMSCF = 1e6/379.49*0.45359237      # 1195.3 кмоль на MMscf (60F, 14.696 psia) [О]
BBL = 0.158987
COMPS = {
 # SPE-184345-MS (Ibemere, Mmata, Onyekonwu, SPE NAICE 2016), Table 2 «Associated Gas Composition», ГХ по ASTM D1945,
 # пробы с flow stations и газовых заводов дельты Нигера; 6 проб [П]. Измеренная HDP проб 18-30 C, GPM C3+ 1.1-2.55 (Fig.2) [П]
 'SPE_AG1': dict(N2=.92,CO2=.54,C1=94.23,C2=3.83,C3=.12,iC4=.13,nC4=.01,iC5=.07,nC5=0,C6=.08,C7=.04,C8=.02,C9=.01,C10=0),
 'SPE_AG2': dict(N2=.03,CO2=1.48,C1=91.08,C2=3.42,C3=2.25,iC4=.42,nC4=.73,iC5=.21,nC5=.12,C6=.11,C7=.08,C8=.05,C9=.02,C10=0),
 'SPE_AG3': dict(N2=.16,CO2=17.44,C1=72.86,C2=4.05,C3=2.58,iC4=.87,nC4=1.00,iC5=.44,nC5=.24,C6=.15,C7=.11,C8=.07,C9=.02,C10=.01),
 'SPE_AG4': dict(N2=.07,CO2=.74,C1=84.07,C2=6.55,C3=5.44,iC4=.93,nC4=1.43,iC5=.26,nC5=.23,C6=.14,C7=.10,C8=.04,C9=0,C10=0),
 'SPE_AG5': dict(N2=.05,CO2=.59,C1=87.83,C2=6.07,C3=3.06,iC4=.54,nC4=.99,iC5=.27,nC5=.21,C6=.17,C7=.13,C8=.09,C9=0,C10=0),
 'SPE_AG6': dict(N2=.03,CO2=.02,C1=90.75,C2=4.28,C3=2.49,iC4=.50,nC4=.89,iC5=.29,nC5=.22,C6=.20,C7=.22,C8=.11,C9=0,C10=0),
 # Marginal Oilfield A, Ndunagu et al., NJTD 19(1) 2022, Table 1 [П] (из lpg_model.py корпуса)
 'NJTD_A':  dict(N2=.10,CO2=.46,C1=82.72,C2=6.73,C3=5.68,iC4=1.56,nC4=1.60,iC5=.40,nC5=.26,C6=.16,C7=.20),
 # тощая проба дельты (IJOGCE 2022) [В], пентаны «малые» -> 0.10/0.08 [Д]
 'lean_ND': dict(C1=87.92,C2=4.65,C3=.93,iC4=.26,nC4=.29,iC5=.10,nC5=.08),
}
c, p = ChemicalConstantsPackage.from_IDs(IDS)
kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', c.CASs, 'kij')
eos = {'Tcs': c.Tcs, 'Pcs': c.Pcs, 'omegas': c.omegas, 'kijs': kijs}
FL = FlashVL(c, p, liquid=CEOSLiquid(PRMIX, eos, HeatCapacityGases=p.HeatCapacityGases),
             gas=CEOSGas(PRMIX, eos, HeatCapacityGases=p.HeatCapacityGases))
def zvec(d):
    v = [d.get(k, 0.0) for k in KEYS]; s = sum(v); return [x/s for x in v]
def liq_bbl(zl, nmol_kmol):
    # объём идеальной жидкости, барр., из кмоль и мольных долей
    m3 = sum(zl[i]*nmol_kmol*c.MWs[i]/(SG[k]*999.0) for i, k in enumerate(KEYS))
    return m3/BBL
def potential(d):
    z = zvec(d); n = KMOL_PER_MMSCF
    c5 = sum(z[i]*n*c.MWs[i]/(SG[k]*999.0) for i, k in enumerate(KEYS) if k in ('iC5','nC5','C6','C7','C8','C9','C10'))/BBL
    c3 = sum(z[i]*n*c.MWs[i]/(SG[k]*999.0) for i, k in enumerate(KEYS) if k in ('C3','iC4','nC4','iC5','nC5','C6','C7','C8','C9','C10'))/BBL
    return dict(C5p_bbl_per_MMscf=round(c5, 2), C3p_bbl_per_MMscf=round(c3, 2),
                GPM_C3p=round(c3*42/1000, 2), GPM_C5p=round(c5*42/1000, 2), C5p_mol=round(100*sum(z[i] for i,k in enumerate(KEYS) if k in ('iC5','nC5','C6','C7','C8','C9','C10')),2))
def sep(z, T_C, P_bar):
    r = FL.flash(T=273.15+T_C, P=P_bar*1e5, zs=z)
    if r.liquid_count == 0 or r.VF >= 1.0: return 0.0, None, r
    return 1.0-r.VF, r.liquid0.zs, r
# Стабилизатор (колонна/подогреваемая ёмкость) моделируется покомпонентным расщеплением в кубовый продукт [Д]:
# C1-C3 уходят в верх (на факел), iC4 15%, nC4 30%, iC5 90%, nC5 95%, C6+ 100% -> ДНП ~0.7-0.8 бар (10-12 psia).
# Одноступенчатый флэш при 45 C / 1.2 бар для жидкости НТС, богатой C3-C4, испаряет и C5+ — это артефакт, не стабилизатор.
SPLIT = dict(C1=0, C2=0, C3=0, iC4=0.15, nC4=0.30, iC5=0.90, nC5=0.95, C6=1, C7=1, C8=1, C9=1, C10=1, N2=0, CO2=0)
def stabilize(xl, L_kmol):
    n = [xl[i]*L_kmol*SPLIT[k] for i, k in enumerate(KEYS)]
    tot = sum(n)
    if tot <= 0: return 0.0, None
    return liq_bbl([v/tot for v in n], tot), [v/tot for v in n]
def cold_sep_yield(d, T_C, P_bar):
    z = zvec(d); L, xl, _ = sep(z, T_C, P_bar)
    if L <= 0: return dict(raw=0.0, stab=0.0, c5_mol_in_stab=0.0)
    Lk = L*KMOL_PER_MMSCF
    raw = liq_bbl(xl, Lk); sb, xs = stabilize(xl, Lk)
    c5frac = 0.0 if xs is None else sum(xs[i] for i,k in enumerate(KEYS) if k in ('iC5','nC5','C6','C7','C8','C9','C10'))
    return dict(raw=round(raw, 2), stab=round(sb, 2), c5_mol_in_stab=round(c5frac, 2))
def dewpoint(d, P_bar):
    z = zvec(d)
    try: return round(FL.flash(P=P_bar*1e5, VF=1.0, zs=z).T-273.15, 1)
    except Exception: return None
def jt_T(d, T1, P1, P2):
    z = zvec(d); h = FL.flash(T=273.15+T1, P=P1*1e5, zs=z).H()
    return round(FL.flash(H=h, P=P2*1e5, zs=z).T-273.15, 1)
def jt_recup(d, T1, P1, P2, pinch=5.0, n=10):
    """Схема: входной сепаратор -> газ/газ теплообменник (противоток) -> холодный сепаратор ВД (P1, Tx) ->
    JT-клапан -> НТС (P2) -> холодный газ НТС греется в теплообменнике -> на факел. Жидкость обоих сепараторов -> стабилизатор.
    Ищем самую низкую Tx (шаг 2 C), при которой минимальный напор в теплообменнике >= pinch.
    Теплопотери и перепад давления в теплообменнике не учтены (оптимистично) [Д]."""
    z = zvec(d); hin = FL.flash(T=273.15+T1, P=P1*1e5, zs=z)
    best = None; Tx = T1 - 2.0
    while Tx > -45:
        hx = FL.flash(T=273.15+Tx, P=P1*1e5, zs=z)
        if hx.liquid_count and hx.gas is not None:
            V1 = hx.VF; y1 = hx.gas.zs; L1 = 1-V1; x1 = hx.liquid0.zs
        else:
            V1 = 1.0; y1 = z; L1 = 0.0; x1 = None
        hv = FL.flash(T=273.15+Tx, P=P1*1e5, zs=y1).H()
        r2 = FL.flash(H=hv, P=P2*1e5, zs=y1)
        Tl = r2.T-273.15
        if r2.liquid_count and r2.gas is not None:
            V2 = V1*r2.VF; y2 = r2.gas.zs; L2 = V1*(1-r2.VF); x2 = r2.liquid0.zs
        else:
            V2 = V1; y2 = y1; L2 = 0.0; x2 = None
        Tg = [Tl + 1.5*i for i in range(int((T1+10-Tl)/1.5)+2)]
        Hg = [V2*FL.flash(T=273.15+t, P=P2*1e5, zs=y2).H() for t in Tg]
        def Tcold(q):
            q = q + Hg[0]
            for i in range(1, len(Hg)):
                if Hg[i] >= q: return Tg[i-1] + 1.5*(q-Hg[i-1])/(Hg[i]-Hg[i-1])
            return Tg[-1]
        ok = True; minDT = 1e9
        for k in range(n+1):
            Th = Tx + (T1-Tx)*k/n
            qh = FL.flash(T=273.15+Th, P=P1*1e5, zs=z).H() - hx.H()
            dT = Th - Tcold(qh); minDT = min(minDT, dT)
            if dT < pinch: ok = False; break
        if not ok: break
        yb = 0.0; parts = []
        for L, x in ((L1, x1), (L2, x2)):
            if L > 0 and x is not None:
                sb, _ = stabilize(x, L*KMOL_PER_MMSCF); yb += sb; parts.append(round(sb, 2))
            else: parts.append(0.0)
        best = dict(Tx=round(Tx, 1), T_LTS=round(Tl, 1), L_HP=round(L1, 5), L_LTS=round(L2, 5), minDT=round(minDT, 1),
                    yield_stab_bbl_per_MMscf=round(yb, 2), split_HP_LTS=parts)
        Tx -= 2.0
    return best
def hydrate_T(d, P_bar):
    # Motiee (1991) по относительной плотности газа [В]; °F <- psia
    z = zvec(d); MW = sum(z[i]*c.MWs[i] for i in range(len(z))); g = MW/28.964
    lp = math.log10(P_bar*14.5038)
    TF = -238.24469 + 78.99667*lp - 5.352544*lp**2 + 349.473877*g - 150.854675*g*g - 27.604065*lp*g
    return round((TF-32)/1.8, 1), round(g, 3)
def water_kg_per_MMscf(T_C, P_bar):
    # насыщение по Раулю (Антуан для воды), грубо +-20% [О]
    psat = 10**(8.07131-1730.63/(233.426+T_C))*133.322/1e5
    return KMOL_PER_MMSCF*psat/P_bar*18.015
def meoh_need(d, T1, P1, T_cold, P_cold):
    Th, _ = hydrate_T(d, P1)
    dT = Th - T_cold + 3.0                                   # запас 3 C [Д]
    if dT <= 0: return dict(hydrate_T=Th, need=False, kg_per_MMscf=0.0)
    W = 100*dT*1.8*32/(2335+dT*1.8*32)                       # Хаммершмидт, мас.% MeOH в водной фазе [В]
    wat = max(0.0, water_kg_per_MMscf(T1, P1) - water_kg_per_MMscf(T_cold, P_cold))
    aq = wat*W/(100-W)
    return dict(hydrate_T=Th, need=True, wt_pct=round(W,1), water_cond=round(wat,1), meoh_aq=round(aq,1),
                kg_per_MMscf_lo=round(aq*1.5,1), kg_per_MMscf_hi=round(aq*4.0,1))   # потери в газ/конденсат x1.5-4 [Д]
if __name__ == '__main__':
    out = dict(potential={}, dew={}, yield_grid={}, jt_adiabatic={}, jt_recup={}, hydrate={})
    for n, d in COMPS.items():
        out['potential'][n] = potential(d)
        out['dew'][n] = {f'{P}bar': dewpoint(d, P) for P in (3, 15, 30, 50)}
        g = {}
        for P in (3, 15, 30, 50):
            for T in (35, 25, 15, 5, -5, -15, -25):
                g[f'{P}bar_{T}C'] = cold_sep_yield(d, T, P)
        out['yield_grid'][n] = g
        out['jt_adiabatic'][n] = {f'{P1}->{P2}bar_from35C': jt_T(d, 35, P1, P2) for P1, P2 in ((15,3),(30,3),(50,3),(50,15),(70,30))}
        out['hydrate'][n] = {f'{P}bar': hydrate_T(d, P)[0] for P in (3, 15, 30, 50)}
    for n in ('SPE_AG2','SPE_AG4','SPE_AG5','SPE_AG6','NJTD_A','lean_ND'):
        rr = {}
        for (P1, P2) in ((15, 3), (30, 3), (50, 3)):
            bst = jt_recup(COMPS[n], 35.0, P1, P2)
            if bst:
                bst['meoh'] = meoh_need(COMPS[n], 35.0, P1, bst['Tx'], P1)
                bst['C5p_recovery_vs_potential'] = round(bst['yield_stab_bbl_per_MMscf']/out['potential'][n]['C5p_bbl_per_MMscf'], 2)
            rr[f'{P1}->{P2}'] = bst
            print(n, P1, P2, bst, flush=True)
        out['jt_recup'][n] = rr
    json.dump(out, open(os.path.join(HERE, 'process_out.json'), 'w'), ensure_ascii=False, indent=1)
    for k in ('potential', 'dew', 'jt_adiabatic', 'hydrate'):
        print('==', k)
        for n, v in out[k].items(): print(n, v)
    print('== yield_grid (stab bbl/MMscf)')
    for n, g in out['yield_grid'].items(): print(n, {k: v['stab'] for k, v in g.items()})
    print('== jt_recup')
    for n, v in out['jt_recup'].items(): print(n, json.dumps(v, ensure_ascii=False))
