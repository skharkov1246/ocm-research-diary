"""Эксергетический анализ двигателя-заменителя KV-3 (трек «Простая литая ГТУ-заменитель», issue #59).

Что считает (уровень облика, 0D):
  * цепочка состояний — точная копия gtu_cycle.cycle_r (Cantera, GRI-Mech 3.0, метан, полное сгорание,
    политропные процессы 40 шагами); совпадение мощности, расхода топлива и выхлопа с cycle/recup
    проверяется в каждом прогоне;
  * эксергия каждого потока = физическая (тепловая + от давления) + химическая; окружающая среда —
    стандартные условия ISO 3977-2 (15 °C, 101,325 кПа, относительная влажность 60 %): сухой воздух модели
    + водяной пар; эксергия метана — через реакцию CH4 + 2 O2 -> CO2 + 2 H2O с химическими потенциалами среды;
  * разрушение эксергии по узлам двумя способами: (A) баланс эксергии узла, (B) теорема Гюи — Стодолы
    T0·S_gen по энтропиям потоков; расхождение A и B — невязка энергии узла;
  * баланс всего двигателя по способу B: эксергия топлива и воздуха = мощность на клеммах + разрушение +
    уходящие потоки (выхлоп, утечки); невязка обязана быть не больше 0,5 % эксергии топлива;
  * три исполнения: KV-3 в расчётной точке (frame3), с рекуператором, с повышенной температурой газа и
    печатным охлаждением; перебор приёмов «где резерв КПД» на расчётной точке;
  * диаграммы Грассмана (поток эксергии) и сравнение узлов — SVG.

Запуск: python3 exergy.py  (около 5 мин на одном ядре) -> exergy_results.json и *.svg рядом со скриптом.
AWS не вызывается. Все КПД узлов, потери и модель охлаждения — оценка (см. assumptions в JSON)."""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO = os.path.abspath(os.path.join(GTU, "..", ".."))
_cwd = os.getcwd()
sys.path.insert(0, GTU)
import gtu_cycle as C  # noqa: E402  (импорт делает chdir в routes/gtu)
import gtu_explore as E  # noqa: E402
import gtu_family as F  # noqa: E402
os.chdir(_cwd)
import cantera as ct  # noqa: E402

FRAME3 = os.path.join(GTU, "frame3", "frame3_point_results.json")
OUT = os.path.join(HERE, "exergy_results.json")
DATE = "2026-11-14"                # дата расчёта (правило трека: date 2026-11-14), не дата перезапуска
T0, P0 = C.T0, C.P0
RH_ISO = 0.60                      # ISO 3977-2: 15 °C, 101,325 кПа, 60 % (по выдержке поиска, см. README)
BALANCE_TOL_PCT = 0.5              # требование задачи: баланс эксергии закрыт до 0,5 % эксергии топлива
MAIN = ("N2", "O2", "AR", "CO2", "H2O")

# ---------------------------------------------------------------- вход: расчётная точка KV-3 (frame3)


def _rel(p):
    return os.path.relpath(p, REPO)


def load_frame3():
    d = json.load(open(FRAME3, encoding="utf-8"))
    inp = {x["field"]: x["used"] for x in d["inputs"]}
    slot = inp["slots[model=Taurus 60, variant=генератор].мягкий"]
    fam = inp["SCEN, COOL['U'], ARCH['1В-2Ц'].eta_pc, ETA_GEN, _run (eta_pt_gg, eta_pt_pt, x_leak)"]
    cyc = d["cycle"]["frame3"]
    t3 = d["turbine"]["options"][str(d["turbine"]["selected_z"])]
    return dict(
        PR=d["split"]["pi"], T41_C=float(slot["T41_C"]), G=float(slot["G"]), n_rpm=float(slot["n_out"]),
        eta_pc=cyc["eta_pc"], eta_t=cyc["eta_t"], dp_duct=cyc["dp_duct"], x_leak=fam["x_leak"],
        eta_gen=fam["eta_gen"], eta_gb=fam["eta_gearbox"],
        P_shaft_kW_ref=cyc["P_shaft_kW"], eff_ref=cyc["eff_shaft"], EGT_ref=cyc["EGT_C"],
        pi1=d["stage1"]["PR"], eta_p1=d["stage1"]["eta_p"], pi2=d["stage2"]["PR"], eta_p2=d["stage2"]["eta_p"],
        eta_p1_nofric=d["stage1"]["eta_basis"]["comp_stage_eta_p"] + min(0.0, d["stage1"]["eta_basis"]["d_eta_without_diffuser_friction"]),
        eta_p2_nofric=d["stage2"]["eta_basis"]["comp_stage_eta_p"] + min(0.0, d["stage2"]["eta_basis"]["d_eta_without_diffuser_friction"]),
        turb_ER=[s["ER"] for s in t3["stages"]], turb_eta=t3["eta_p"],
        T_rel1_C=t3["stages"][0]["T_rel_rotor_C"], T01_1_C=t3["stages"][0]["T01_C"],
        P_req_kW=d["cycle"]["required"]["P_shaft_kW"], G_cap=d["cycle"]["required"]["G_cap"],
        EGT_cap=d["cycle"]["required"]["EGT_cap"],
    )


# ---------------------------------------------------------------- окружающая среда и эксергия потока
class Env:
    """Среда ISO: сухой воздух модели (gtu_cycle.AIR) + водяной пар при RH; химические потенциалы среды."""

    def __init__(self, T=T0, P=P0, RH=RH_ISO, dry=C.AIR, x_override=None):
        self.T, self.P = T, P
        g = ct.Solution("gri30.yaml")
        self.g = g
        self.M = g.molecular_weights
        self.idx = {k: g.species_index(k) for k in MAIN + ("CH4",)}
        if x_override is None:
            w = ct.Water()
            w.TQ = T, 0.0
            self.p_sat = w.P
            self.x_w = RH * self.p_sat / P
            g.TPX = T, P, dry
            X = g.X * (1.0 - self.x_w)
            X[self.idx["H2O"]] += self.x_w
        else:
            self.p_sat, self.x_w = None, x_override.get("H2O", 0.0)
            X = np.zeros(g.n_species)
            for k, v in x_override.items():
                X[self.idx[k]] = v
            X /= X.sum()
        g.TPX = T, P, X
        self.X = {k: float(g.X[self.idx[k]]) for k in MAIN}
        mu = g.chemical_potentials                       # Дж/кмоль
        self.mu00 = {k: float(mu[self.idx[k]]) for k in MAIN}
        self.mu00["CH4"] = self.mu00["CO2"] + 2 * self.mu00["H2O"] - 2 * self.mu00["O2"]
        self.vec = np.full(g.n_species, np.nan)
        for k, v in self.mu00.items():
            self.vec[self.idx[k]] = v

    def ex(self, T, P, Y):
        """Эксергия, Дж/кг: тепловая, от давления, химическая (состав заморожен; кинетическая не учитывается)."""
        g, T0_, P0_ = self.g, self.T, self.P
        Y = np.asarray(Y)
        g.TPY = T, P, Y
        h, s = g.enthalpy_mass, g.entropy_mass
        g.TPY = T0_, P, Y
        hp, sp = g.enthalpy_mass, g.entropy_mass
        g.TPY = T0_, P0_, Y
        h0, s0 = g.enthalpy_mass, g.entropy_mass
        mu = g.chemical_potentials
        nz = Y > 1e-14
        if np.isnan(self.vec[nz]).any():
            bad = [g.species_name(i) for i in np.where(nz & np.isnan(self.vec))[0]]
            raise ValueError(f"нет химического потенциала среды для {bad}")
        ch = float(np.sum(Y[nz] / self.M[nz] * (mu[nz] - self.vec[nz])))
        exT = (h - hp) - T0_ * (s - sp)
        exP = (hp - h0) - T0_ * (sp - s0)
        return dict(T=exT, P=exP, ch=ch, tot=exT + exP + ch, h=h, s=s)


def fuel_exergy_check(env):
    """Эксергия метана двумя способами (одни и те же данные NASA-полиномов gri30, разная алгебра):
    (1) химические потенциалы среды (как в расчёте); (2) -ΔH_r + T0·ΔS_r реакции чистых веществ при T0, p0
    плюс RT0·ln(p0/p_i) для O2, CO2, H2O среды. Плюс расчёт в среде Шаргута (298,15 К) для сверки с таблицей."""
    g = ct.Solution("gri30.yaml")
    R = ct.gas_constant
    hs = {}
    for k in ("CH4", "O2", "CO2", "H2O"):
        g.TPX = env.T, env.P, f"{k}:1"
        hs[k] = (g.enthalpy_mole, g.entropy_mole)
    dH = hs["CO2"][0] + 2 * hs["H2O"][0] - hs["CH4"][0] - 2 * hs["O2"][0]
    dS = hs["CO2"][1] + 2 * hs["H2O"][1] - hs["CH4"][1] - 2 * hs["O2"][1]
    m2 = -(dH - env.T * dS) + R * env.T * (math.log(1 / env.X["CO2"]) + 2 * math.log(1 / env.X["H2O"]) - 2 * math.log(1 / env.X["O2"]))
    g.TPX = env.T, env.P, "CH4:1"
    m1 = g.chemical_potentials[env.idx["CH4"]] - env.mu00["CH4"]
    # теплота сгорания из тех же данных (продукты - пар, при T0): низшая
    lhv = -dH / env.M[env.idx["CH4"]]
    sz = Env(T=298.15, P=101325.0, x_override=dict(N2=75.78, O2=20.39, CO2=0.0335, H2O=2.2, AR=0.906))
    g.TPX = 298.15, 101325.0, "CH4:1"
    m3 = g.chemical_potentials[sz.idx["CH4"]] - sz.mu00["CH4"]
    return dict(ex_ch_kJ_mol_mu=m1 / 1e6, ex_ch_kJ_mol_reaction=m2 / 1e6, diff_rel=abs(m1 - m2) / m1,
                ex_ch_MJ_kg=m1 / env.M[env.idx["CH4"]] / 1e6, LHV_cantera_MJ_kg=lhv / 1e6, LHV_model_MJ_kg=C.LHV / 1e6,
                ratio_ex_LHV_model=m1 / env.M[env.idx["CH4"]] / C.LHV,
                szargut_env_298K_kJ_mol=m3 / 1e6,
                szargut_env_note="среда Шаргута (парциальные давления, кПа: N2 75,78; O2 20,39; CO2 0,0335; H2O 2,2; Ar 0,906) — по памяти, источник не открыт; "
                                 "табличное значение 831,65 кДж/моль (Szargut, Morris, Steward 1988) выдержкой поиска не подтверждено")


# ---------------------------------------------------------------- цепочка состояний (копия gtu_cycle.cycle_r)
def _snap(g):
    return dict(T=g.T, P=g.P, Y=g.Y.copy(), h=g.enthalpy_mass, s=g.entropy_mass)


def _state(T=None, P=None, Y=None, H=None):
    g = ct.Solution("gri30.yaml")
    if H is None:
        g.TPY = T, P, Y
    else:
        g.HPY = H, P, Y
    return _snap(g)


def chain(PR, T41_C, eta_pc, eta_pt_gg, eta_pt_pt, x_ngv, x_rot, x_leak, dp_comb=0.04, dp_duct=0.01,
          p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, W1=1.0, T3r=None, dp_rec_air=0.0, T_fuel=None):
    """Та же последовательность операций и те же допуски brentq, что в gtu_cycle.cycle_r; дополнительно
    сохраняет все промежуточные состояния. T_fuel — расширение (подогрев топлива); None = T0, как в модели."""
    T41 = T41_C + 273.15
    Tf = T0 if T_fuel is None else T_fuel
    a = C.gas_air()
    st1 = _snap(a)
    wc = C.poly_process(a, P0 * PR, eta_pc, True)
    st3 = _snap(a)
    T3, p3 = a.T, a.P
    Wca = W1 * (1 - x_ngv - x_rot - x_leak)

    def after_mix(far, keep=False):
        g = ct.Solution("gri30.yaml")
        f = ct.Quantity(g, constant="HP"); f.TPX = Tf, p3, "CH4:1"; f.mass = far * Wca
        ai = ct.Quantity(g, constant="HP"); ai.TPX = (T3r or T3), p3, C.AIR; ai.mass = Wca
        sf = dict(T=f.T, P=f.P, Y=f.Y.copy(), h=f.enthalpy_mass, s=f.entropy_mass)
        sa = dict(T=ai.T, P=ai.P, Y=ai.Y.copy(), h=ai.enthalpy_mass, s=ai.entropy_mass)
        m = ai + f
        hm = (ai.enthalpy_mass * ai.mass + f.enthalpy_mass * f.mass) / (ai.mass + f.mass)
        Ym = m.Y.copy()
        Y = m.Y.copy(); gi = g.species_index
        MW = g.molecular_weights
        nf = Y[gi("CH4")] / MW[gi("CH4")]
        Y[gi("CH4")] = 0
        Y[gi("O2")] -= 2 * nf * MW[gi("O2")]
        Y[gi("CO2")] += nf * MW[gi("CO2")]
        Y[gi("H2O")] += 2 * nf * MW[gi("H2O")]
        p4 = p3 * (1 - dp_comb) * (1 - dp_rec_air)
        Q = (1 - eta_b) * far / (1 + far) * C.LHV
        g.HPY = hm - Q, p4, Y
        T4 = g.T
        st4 = _snap(g)
        q4 = ct.Quantity(g, constant="HP"); q4.mass = Wca * (1 + far)
        g2 = ct.Solution("gri30.yaml")
        c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p4, C.AIR; c.mass = W1 * x_ngv
        mix = q4 + c
        if not keep:
            return mix, T4
        return mix, T4, dict(fuel=sf, air=sa, mix_unreacted=_state(P=p3, Y=Ym, H=hm),
                             prod_ad_p3=_state(P=p3, Y=Y, H=hm), prod_ad_p4=_state(P=p4, Y=Y, H=hm),
                             st4=st4, cool_ngv=_state(T=T3, P=p4, Y=c.Y.copy()), q_hl=Q * Wca * (1 + far), p4=p4)

    far = brentq(lambda f: after_mix(f)[0].T - T41, 0.005, 0.05, xtol=1e-7)
    mix, T4, comb = after_mix(far, keep=True)
    g = ct.Solution("gri30.yaml"); g.TPY = mix.T, mix.P, mix.Y
    st41 = _snap(g)
    W41 = mix.mass
    Wf = far * Wca
    need = wc * W1 / eta_m / W41

    def gg_resid(p_out):
        g.TPY = mix.T, mix.P, mix.Y
        return C.poly_process(g, p_out, eta_pt_gg, False) - need
    p45 = brentq(gg_resid, P0 * 1.05, mix.P * 0.95, xtol=1.0)
    w_gg = gg_resid(p45) + need
    st45 = _snap(g)
    q = ct.Quantity(g, constant="HP"); q.mass = W41
    g2 = ct.Solution("gri30.yaml")
    c = ct.Quantity(g2, constant="HP"); c.TPX = T3, p45, C.AIR; c.mass = W1 * x_rot
    st_cr = dict(T=c.T, P=c.P, Y=c.Y.copy(), h=c.enthalpy_mass, s=c.entropy_mass)
    m2 = q + c
    st45m = _state(T=m2.T, P=m2.P, Y=m2.Y)
    g.TPY = m2.T, p45 * (1 - dp_duct), m2.Y
    st45d = _snap(g)
    W45 = m2.mass
    w_pt = C.poly_process(g, P0 * p_exh_ratio, eta_pt_pt, False)
    st5 = _snap(g)
    P_shaft = w_pt * W45 * eta_m
    return dict(
        inp=dict(PR=PR, T41_C=T41_C, eta_pc=eta_pc, eta_pt_gg=eta_pt_gg, eta_pt_pt=eta_pt_pt, x_ngv=x_ngv, x_rot=x_rot,
                 x_leak=x_leak, dp_comb=dp_comb, dp_duct=dp_duct, p_exh_ratio=p_exh_ratio, eta_b=eta_b, eta_m=eta_m,
                 W1=W1, T3r=T3r, dp_rec_air=dp_rec_air, T_fuel=Tf),
        m=dict(W1=W1, Wca=Wca, Wf=Wf, M=Wca + Wf, W41=W41, W45=W45, Wngv=W1 * x_ngv, Wrot=W1 * x_rot, Wleak=W1 * x_leak),
        st=dict(s1=st1, s3=st3, fuel=comb["fuel"], air_in=comb["air"], mix=comb["mix_unreacted"], pa=comb["prod_ad_p3"],
                pb=comb["prod_ad_p4"], s4=comb["st4"], cngv=comb["cool_ngv"], s41=st41, s45=st45, crot=st_cr,
                s45m=st45m, s45d=st45d, s5=st5),
        w=dict(wc=wc, W_comp=wc * W1, W_gg=w_gg * W41, W_pt=w_pt * W45, P_shaft=P_shaft, Q_hl=comb["q_hl"]),
        out=dict(P_MW=P_shaft / 1e6, eff=P_shaft / (Wf * C.LHV), EGT_C=st5["T"] - 273.15, T3_C=T3 - 273.15,
                 T4_C=T4 - 273.15, T45_C=st45d["T"] - 273.15, Wf=Wf, W_exh=W45, p45=p45, far=far),
    )


def compare_with_model(ch, ref):
    """Цепочка обязана совпасть с gtu_cycle.cycle/recup: относительные расхождения."""
    o = ch["out"]
    return dict(P=abs(o["P_MW"] - ref["P_MW"]) / ref["P_MW"], Wf=abs(o["Wf"] - ref["Wf"]) / ref["Wf"],
                EGT_K=abs(o["EGT_C"] - ref["EGT_C"]), eff=abs(o["eff"] - ref["eff"]))


# ---------------------------------------------------------------- учёт эксергии
GROUPS = [  # ключ группы, подпись, узлы
    ("compressor", "компрессор", ("compressor",)),
    ("comb_mixing", "камера: смешение", ("comb_mixing",)),
    ("comb_reaction", "камера: горение", ("comb_reaction",)),
    ("comb_dp", "камера: потеря давления", ("comb_dp",)),
    ("comb_heat_loss", "камера: недожог", ("comb_heat_loss",)),
    ("turbine", "турбина", ("turbine_gg", "turbine_pt", "duct_dp")),
    ("leak_cool", "утечки и охлаждение", ("leakage", "cool_throttle_ngv", "cool_mix_ngv", "cool_throttle_rotor", "cool_mix_rotor")),
    ("recuperator", "рекуператор", ("recup_heat", "recup_dp_air", "recup_dp_gas")),
    ("exhaust", "выхлоп", ("exhaust",)),
    ("mech_gen", "механика, редуктор, генератор", ("mech", "gearbox", "generator")),
]
LABEL = {
    "compressor": "компрессор (обе ступени)", "comb_mixing": "камера: смешение топлива с воздухом",
    "comb_reaction": "камера: горение (реакция)", "comb_dp": "камера: потеря давления",
    "comb_heat_loss": "камера: недожог η_b (потеря тепла в модели)", "cool_mix_ngv": "охлаждение соплового аппарата: смешение",
    "cool_throttle_ngv": "охлаждение соплового аппарата: дросселирование воздуха p3 -> p4",
    "cool_throttle_rotor": "охлаждение ротора: дросселирование воздуха p3 -> p45",
    "turbine_gg": "турбина: часть, вращающая компрессор", "cool_mix_rotor": "охлаждение ротора: смешение",
    "duct_dp": "переходный канал: потеря давления", "turbine_pt": "турбина: часть, дающая мощность",
    "leakage": "утечки (сжатый воздух уходит мимо турбины)", "recup_heat": "рекуператор: теплообмен при конечной разности температур",
    "recup_dp_air": "рекуператор: потеря давления воздуха", "recup_dp_gas": "рекуператор: потеря давления газа",
    "exhaust": "выхлоп (уходит с газом)", "mech": "механические потери (подшипники)", "gearbox": "редуктор пакета",
    "generator": "генератор пакета",
}
LOSS_NODES = ("leakage", "exhaust")   # эксергия уходит из двигателя с потоком (разрушается вне его)


def account(ch, env, eta_gb, eta_gen, recup=None, split=None):
    """Разрушение эксергии по узлам: способ A (баланс эксергии узла) и B (T0·S_gen). Возвращает словарь узлов,
    балансы и доли от эксергии топлива. recup = dict(dp_gas=..., dp_air=...) для исполнения с рекуператором."""
    T0_ = env.T
    st, m, w = ch["st"], ch["m"], ch["w"]
    X = {k: env.ex(v["T"], v["P"], v["Y"]) for k, v in st.items()}
    Ex = lambda k, mass: X[k]["tot"] * mass  # noqa: E731
    S = lambda k, mass: st[k]["s"] * mass    # noqa: E731
    H = lambda k, mass: st[k]["h"] * mass    # noqa: E731
    nodes = {}

    def put(key, A, B, dH=0.0, S_gen=None):
        nodes[key] = dict(A_W=A, B_W=B, S_gen_W_K=(B / T0_ if S_gen is None else S_gen), dE_W=dH)

    eta_m = ch["inp"]["eta_m"]
    W1, Wca, Wf, M, W41, W45 = m["W1"], m["Wca"], m["Wf"], m["M"], m["W41"], m["W45"]
    # компрессор
    put("compressor", w["W_comp"] + Ex("s1", W1) - Ex("s3", W1), T0_ * (S("s3", W1) - S("s1", W1)),
        dH=w["W_comp"] - (H("s3", W1) - H("s1", W1)))
    # утечка: сжатый воздух уходит из двигателя (его эксергия — потеря)
    nodes["leakage"] = dict(A_W=Ex("s3", m["Wleak"]), B_W=Ex("s3", m["Wleak"]), S_gen_W_K=0.0, dE_W=0.0, outflow=True)
    # рекуператор: холодная сторона s3 -> air_in (T3r, p3); горячая s5 -> стек (баланс энергии)
    st6 = None
    if recup is not None:
        h6 = st["s5"]["h"] - Wca * (st["air_in"]["h"] - st["s3"]["h"]) / W45
        p6 = st["s5"]["P"] / (1 + recup["dp_gas"])
        st6 = _state(P=p6, Y=st["s5"]["Y"], H=h6)
        st["s6"] = st6
        X["s6"] = env.ex(st6["T"], st6["P"], st6["Y"])
        s5p = _state(P=p6, Y=st["s5"]["Y"], H=st["s5"]["h"])    # газ после потери давления при той же энтальпии
        X["s5p"] = env.ex(s5p["T"], s5p["P"], s5p["Y"])
        st["s5p"] = s5p
        A_dp_gas = W45 * (X["s5"]["tot"] - X["s5p"]["tot"])
        B_dp_gas = T0_ * W45 * (s5p["s"] - st["s5"]["s"])
        A_ht = Wca * (X["s3"]["tot"] - X["air_in"]["tot"]) + W45 * (X["s5p"]["tot"] - X["s6"]["tot"])
        B_ht = T0_ * (Wca * (st["air_in"]["s"] - st["s3"]["s"]) + W45 * (st6["s"] - s5p["s"]))
        dHr = Wca * (st["air_in"]["h"] - st["s3"]["h"]) - W45 * (st["s5"]["h"] - h6)
        put("recup_heat", A_ht, B_ht, dH=dHr)
        put("recup_dp_gas", A_dp_gas, B_dp_gas)
    # камера: смешение (воздух + топливо, адиабатно, p3), горение (адиабатно, p3), потеря давления, недожог
    E_air_in = Ex("air_in", Wca)
    put("comb_mixing", E_air_in + Ex("fuel", Wf) - Ex("mix", M),
        T0_ * (S("mix", M) - S("air_in", Wca) - S("fuel", Wf)), dH=H("air_in", Wca) + H("fuel", Wf) - H("mix", M))
    put("comb_reaction", M * (X["mix"]["tot"] - X["pa"]["tot"]), T0_ * M * (st["pa"]["s"] - st["mix"]["s"]),
        dH=M * (st["mix"]["h"] - st["pa"]["h"]))
    A_dp = M * (X["pa"]["tot"] - X["pb"]["tot"])
    B_dp = T0_ * M * (st["pb"]["s"] - st["pa"]["s"])
    dp_c, dp_ra = ch["inp"]["dp_comb"], ch["inp"]["dp_rec_air"]
    if dp_ra > 0:   # в модели потеря давления воздуха рекуператора отнесена на выход камеры: делим по ln
        k = math.log(1 / (1 - dp_ra)) / (math.log(1 / (1 - dp_ra)) + math.log(1 / (1 - dp_c)))
        put("recup_dp_air", A_dp * k, B_dp * k)
        A_dp, B_dp = A_dp * (1 - k), B_dp * (1 - k)
    put("comb_dp", A_dp, B_dp)
    Q = w["Q_hl"]   # теплота, отнятая моделью у продуктов (недожог η_b); граница — среда при T0
    put("comb_heat_loss", M * (X["pb"]["tot"] - X["s4"]["tot"]), T0_ * M * (st["s4"]["s"] - st["pb"]["s"]) + Q,
        dH=M * (st["pb"]["h"] - st["s4"]["h"]) - Q)
    # охлаждающий воздух: дросселирование от p3 до давления подмешивания (в модели — неявно, h = const)
    put("cool_throttle_ngv", Ex("s3", m["Wngv"]) - Ex("cngv", m["Wngv"]), T0_ * m["Wngv"] * (st["cngv"]["s"] - st["s3"]["s"]),
        dH=m["Wngv"] * (st["s3"]["h"] - st["cngv"]["h"]))
    put("cool_throttle_rotor", Ex("s3", m["Wrot"]) - Ex("crot", m["Wrot"]), T0_ * m["Wrot"] * (st["crot"]["s"] - st["s3"]["s"]),
        dH=m["Wrot"] * (st["s3"]["h"] - st["crot"]["h"]))
    # охлаждение соплового аппарата: смешение при p4
    put("cool_mix_ngv", Ex("s4", M) + Ex("cngv", m["Wngv"]) - Ex("s41", W41),
        T0_ * (S("s41", W41) - S("s4", M) - S("cngv", m["Wngv"])), dH=H("s4", M) + H("cngv", m["Wngv"]) - H("s41", W41))
    # турбина (часть, вращающая компрессор)
    put("turbine_gg", Ex("s41", W41) - Ex("s45", W41) - w["W_gg"], T0_ * W41 * (st["s45"]["s"] - st["s41"]["s"]),
        dH=W41 * (st["s41"]["h"] - st["s45"]["h"]) - w["W_gg"])
    put("cool_mix_rotor", Ex("s45", W41) + Ex("crot", m["Wrot"]) - Ex("s45m", W45),
        T0_ * (S("s45m", W45) - S("s45", W41) - S("crot", m["Wrot"])), dH=H("s45", W41) + H("crot", m["Wrot"]) - H("s45m", W45))
    put("duct_dp", W45 * (X["s45m"]["tot"] - X["s45d"]["tot"]), T0_ * W45 * (st["s45d"]["s"] - st["s45m"]["s"]),
        dH=W45 * (st["s45m"]["h"] - st["s45d"]["h"]))
    put("turbine_pt", Ex("s45d", W45) - Ex("s5", W45) - w["W_pt"], T0_ * W45 * (st["s5"]["s"] - st["s45d"]["s"]),
        dH=W45 * (st["s45d"]["h"] - st["s5"]["h"]) - w["W_pt"])
    # механика (на обеих частях турбины, как в модели), редуктор и генератор пакета
    mech = (w["W_gg"] + w["W_pt"]) * (1 - eta_m)
    gb = w["P_shaft"] * (1 - eta_gb)
    gen = w["P_shaft"] * eta_gb * (1 - eta_gen)
    P_el = w["P_shaft"] * eta_gb * eta_gen
    for k_, v in (("mech", mech), ("gearbox", gb), ("generator", gen)):
        put(k_, v, v)
    # невязка мощности вала газогенератора: работа турбины·η_m против работы компрессора (допуск brentq)
    shaft_resid = w["W_gg"] * eta_m - w["W_comp"]
    # выхлоп
    sx = "s6" if st6 is not None else "s5"
    Xe = X[sx]
    nodes["exhaust"] = dict(A_W=Xe["tot"] * W45, B_W=Xe["tot"] * W45, S_gen_W_K=0.0, dE_W=0.0, outflow=True,
                            thermal_W=Xe["T"] * W45, pressure_W=Xe["P"] * W45, chemical_W=Xe["ch"] * W45)
    # баланс
    E_fuel = Ex("fuel", Wf)
    E_air = Ex("s1", W1)
    E_in = E_fuel + E_air
    sumA = sum(v["A_W"] for v in nodes.values())
    sumB = sum(v["B_W"] for v in nodes.values())
    resA = E_in - P_el - sumA
    resB = E_in - P_el - sumB
    # баланс энергии всего двигателя (энтальпия; топливо и воздух — с теплотой образования Cantera)
    H_in = H("fuel", Wf) + H("s1", W1)
    H_out = st[sx]["h"] * W45 + st["s3"]["h"] * m["Wleak"] + P_el + mech + gb + gen + Q
    dE = H_in - H_out
    destr = {k: v for k, v in nodes.items() if not v.get("outflow")}
    D_total = sum(v["B_W"] for v in destr.values())
    S_gen_total = sum(v["S_gen_W_K"] for v in destr.values())
    pct = lambda x: 100.0 * x / E_fuel  # noqa: E731
    groups = []
    for gk, gl, keys in GROUPS:
        val = sum(nodes[k]["B_W"] for k in keys if k in nodes)
        sg = sum(nodes[k]["S_gen_W_K"] for k in keys if k in nodes)
        groups.append(dict(key=gk, label=gl, kW=val / 1e3, pct_fuel_ex=pct(val), S_gen_kW_K=sg / 1e3,
                           kind="уходит с потоком" if gk == "exhaust" else ("разрушение + утечки" if gk == "leak_cool" else "разрушение")))
    res = dict(
        fuel=dict(Ex_kW=E_fuel / 1e3, ex_MJ_kg=X["fuel"]["tot"] / 1e6, ex_ch_MJ_kg=X["fuel"]["ch"] / 1e6,
                  ex_pressure_MJ_kg=X["fuel"]["P"] / 1e6, ex_thermal_MJ_kg=X["fuel"]["T"] / 1e6,
                  LHV_kW=Wf * C.LHV / 1e3, ratio_ex_LHV=X["fuel"]["tot"] / C.LHV),
        air_in=dict(Ex_kW=E_air / 1e3, ex_kJ_kg=X["s1"]["tot"] / 1e3,
                    note="сухой воздух модели относительно влажной среды ISO: химическая эксергия"),
        P_shaft_kW=w["P_shaft"] / 1e3, P_el_kW=P_el / 1e3,
        eff_shaft=w["P_shaft"] / (Wf * C.LHV), eff_el=P_el / (Wf * C.LHV),
        ex_eff_shaft=w["P_shaft"] / E_fuel, ex_eff_el=P_el / E_fuel,
        nodes={k: dict(label=LABEL[k], kW=v["B_W"] / 1e3, kW_by_exergy_balance=v["A_W"] / 1e3, pct_fuel_ex=pct(v["B_W"]),
                       S_gen_kW_K=v["S_gen_W_K"] / 1e3, A_minus_B_pct_fuel_ex=pct(v["A_W"] - v["B_W"]),
                       energy_residual_kW=v["dE_W"] / 1e3,
                       **({"outflow": True} if v.get("outflow") else {}),
                       **({k2: v[k2] / 1e3 for k2 in ("thermal_W", "pressure_W", "chemical_W") if k2 in v}))
               for k, v in nodes.items()},
        groups=groups,
        totals=dict(destruction_kW=D_total / 1e3, destruction_pct=pct(D_total), S_gen_kW_K=S_gen_total / 1e3,
                    entropy_generation_number=D_total / E_fuel,
                    outflow_kW=(nodes["exhaust"]["B_W"] + nodes["leakage"]["B_W"]) / 1e3,
                    outflow_pct=pct(nodes["exhaust"]["B_W"] + nodes["leakage"]["B_W"]), product_pct=pct(P_el),
                    air_in_pct=pct(E_air)),
        balance=dict(E_in_kW=E_in / 1e3, residual_B_kW=resB / 1e3, residual_B_pct_fuel_ex=pct(resB),
                     residual_A_pct_fuel_ex=pct(resA), energy_residual_kW=dE / 1e3,
                     energy_residual_pct_LHV=100 * dE / (Wf * C.LHV), gg_shaft_residual_kW=shaft_resid / 1e3,
                     closed=abs(pct(resB)) <= BALANCE_TOL_PCT,
                     method="B: эксергия топлива + воздуха = мощность на клеммах + Σ T0·S_gen узлов + эксергия выхлопа и утечек; "
                            "A — то же с разрушением из баланса эксергии каждого узла (тождество, контроль арифметики)"),
        states={k: dict(T_C=v["T"] - 273.15, p_kPa=v["P"] / 1e3, ex_kJ_kg=X[k]["tot"] / 1e3, ex_ch_kJ_kg=X[k]["ch"] / 1e3)
                for k, v in st.items()},
        flows=dict(W1=m["W1"], Wca=Wca, Wf=Wf, W41=W41, W45=W45, Wngv=m["Wngv"], Wrot=m["Wrot"], Wleak=m["Wleak"]),
        works_kW=dict(W_comp=w["W_comp"] / 1e3, W_turb=(w["W_gg"] + w["W_pt"]) / 1e3, W_gg=w["W_gg"] / 1e3, W_pt=w["W_pt"] / 1e3),
        Ex_kW=dict(air_c3=Ex("s3", Wca) / 1e3, air_comb_in=E_air_in / 1e3, cool=Ex("s3", m["Wngv"] + m["Wrot"]) / 1e3,
                   fuel=E_fuel / 1e3, s4=Ex("s4", M) / 1e3, s41=Ex("s41", W41) / 1e3, s5=Ex("s5", W45) / 1e3,
                   exhaust=Xe["tot"] * W45 / 1e3, recup_transfer=(E_air_in - Ex("s3", Wca)) / 1e3 if recup else 0.0),
    )
    if split is not None:
        res["compressor_stages"] = split
    if st6 is not None:
        Th, Tc = st["s5"]["T"], st["s3"]["T"]
        q_kW = Wca * (st["air_in"]["h"] - st["s3"]["h"]) / 1e3
        dT1, dT2 = Th - st["air_in"]["T"], st6["T"] - Tc
        lmtd = (dT1 - dT2) / math.log(dT1 / dT2) if abs(dT1 - dT2) > 1e-6 else dT1
        res["recuperator"] = dict(Q_kW=q_kW, T_gas_in_C=Th - 273.15, T_gas_out_C=st6["T"] - 273.15,
                                  T_air_in_C=Tc - 273.15, T_air_out_C=st["air_in"]["T"] - 273.15,
                                  eps_actual=(st["air_in"]["T"] - Tc) / (Th - Tc), LMTD_K=lmtd, UA_kW_K=q_kW / lmtd)
    return res


def compressor_split(ch, pi1, eta_p1):
    """Разбивка компрессора на две ступени: ступень 1 — π и КПД frame3, ступень 2 — КПД, при котором выход
    совпадает с выходом цепочки (та же энтальпия при той же π): разрушение ступеней T0·W1·Δs."""
    st1, st3, W1 = ch["st"]["s1"], ch["st"]["s3"], ch["m"]["W1"]
    if ch["inp"]["PR"] <= pi1 * (1 + 1e-9):
        return dict(n_stages=1, note="π не больше π ступени 1 frame3: одна ступень")
    a = C.gas_air()
    C.poly_process(a, P0 * pi1, eta_p1, True)
    s2 = _snap(a)

    def resid(eta):
        g = ct.Solution("gri30.yaml"); g.TPY = s2["T"], s2["P"], s2["Y"]
        C.poly_process(g, st3["P"], eta, True)
        return g.enthalpy_mass - st3["h"]
    eta2 = brentq(resid, 0.5, 0.999, xtol=1e-9)
    D1 = T0 * W1 * (s2["s"] - st1["s"])
    D2 = T0 * W1 * (st3["s"] - s2["s"])
    return dict(n_stages=2, pi1=pi1, eta_p1=eta_p1, pi2=ch["inp"]["PR"] / pi1, eta_p2_matched=eta2, T2_C=s2["T"] - 273.15,
                D1_kW=D1 / 1e3, D2_kW=D2 / 1e3)


def turbine_stages(ch, ERs, eta):
    """Разрушение по ступеням турбины (frame3: три ступени, степени расширения пропорционально в логарифме
    приведены к степени расширения цепочки). Только для схемы без охлаждения ротора и без канала."""
    s41, s5, W = ch["st"]["s41"], ch["st"]["s5"], ch["m"]["W41"]
    ER_tot = s41["P"] / s5["P"]
    k = math.log(ER_tot) / sum(math.log(x) for x in ERs)
    g = ct.Solution("gri30.yaml"); g.TPY = s41["T"], s41["P"], s41["Y"]
    out = []
    for i, er in enumerate(ERs):
        sa = g.entropy_mass; Ta = g.T
        w = C.poly_process(g, g.P / er ** k, eta, False)
        out.append(dict(stage=i + 1, ER=er ** k, T_in_C=Ta - 273.15, T_out_C=g.T - 273.15, w_kJkg=w / 1e3,
                        D_kW=T0 * W * (g.entropy_mass - sa) / 1e3))
    return dict(stages=out, sum_D_kW=sum(x["D_kW"] for x in out), ER_scale_log=k)


# ---------------------------------------------------------------- модель охлаждения (печатное против литого)
def cooling_model(f3, T3_C):
    """Доля охлаждающего воздуха ряда: x = K·ε/(1-ε), ε = (T_газа - T_стенки)/(T_газа - T_охл) — структура
    Holland & Thake (1980), Horlock et al. (2001); K ∝ 1/η_c (η_c — эффективность внутреннего охлаждения).
    Калибровка литья: класс K трека (gtu_family.COOL['K']: 5,5 % + 5,5 % при T41 1200 °C) при T3 расчётной π.
    Допустимая температура: сопловой аппарат — T41 класса U (газ до 1000 °C без охлаждения), ротор — та же граница
    в относительной системе (frame3: ступень 1, T_rel = T01 - dT_rel). Печать: η_c 0,75 против 0,5 у литья (оценка)."""
    U = F.COOL["U"]; Kc = F.COOL["K"]
    Tm_ngv = float(max(U["T"]))
    dT_rel = f3["T01_1_C"] - f3["T_rel1_C"]
    Tm_rot = Tm_ngv - dT_rel
    Tcal = float(max(Kc["T"]))

    def phi(Tg, Tm, Tc):
        e = (Tg - Tm) / (Tg - Tc)
        return 0.0 if e <= 0 else e / (1 - e)
    K_ngv = Kc["xn"] / phi(Tcal, Tm_ngv, T3_C)
    K_rot = Kc["xr"] / phi(Tcal - dT_rel, Tm_rot, T3_C)
    eta_c = dict(cast=0.5, printed=0.75)

    def x(T41_C, tech, T3=T3_C):
        f = eta_c["cast"] / eta_c[tech]
        return (f * K_ngv * phi(T41_C, Tm_ngv, T3), f * K_rot * phi(T41_C - dT_rel, Tm_rot, T3))
    return x, dict(T_m_ngv_C=Tm_ngv, dT_rel_K=dT_rel, T_m_rot_C=Tm_rot, T41_cal_C=Tcal, x_cal=(Kc["xn"], Kc["xr"]),
                   T3_cal_C=T3_C, K_ngv_cast=K_ngv, K_rot_cast=K_rot, eta_c=eta_c,
                   track_classes={k: dict(T_max=max(v["T"]), xn=v["xn"], xr=v["xr"]) for k, v in F.COOL.items()})


# ---------------------------------------------------------------- исполнения
def base_kw(f3, W1=None):
    return dict(PR=f3["PR"], T41_C=f3["T41_C"], eta_pc=f3["eta_pc"], eta_pt_gg=f3["eta_t"], eta_pt_pt=f3["eta_t"],
                x_ngv=0.0, x_rot=0.0, x_leak=f3["x_leak"], dp_duct=f3["dp_duct"], W1=f3["G"] if W1 is None else W1)


def _scaled_ref(ref, kw, P_target_kW):
    """Модель линейна по расходу: при P_target пересчитываем W1 и экстенсивные величины эталона."""
    if P_target_kW is None:
        return ref, kw
    k = P_target_kW / (ref["P_MW"] * 1e3)
    ref = dict(ref, P_MW=ref["P_MW"] * k, Wf=ref["Wf"] * k, W_exh=ref["W_exh"] * k)
    return ref, dict(kw, W1=kw["W1"] * k)


def run_simple(kw, env, f3, split=True, P_target=None):
    ref = C.cycle(**{k: v for k, v in kw.items() if k != "T_fuel"}) if kw.get("T_fuel") is None else None
    if ref is not None:
        ref, kw = _scaled_ref(ref, kw, P_target)
    ch = chain(**kw)
    acc = account(ch, env, f3["eta_gb"], f3["eta_gen"],
                  split=compressor_split(ch, f3["pi1"], f3["eta_p1"]) if split else None)
    acc["model_check"] = compare_with_model(ch, ref) if ref else None
    acc["inputs"] = {k: v for k, v in kw.items()}
    acc["cycle"] = {k: v for k, v in ch["out"].items()}
    return acc, ch


def run_recup(kw, env, f3, eps, dp_air=0.025, dp_gas=0.035, P_target=None):
    r = C.recup(eps=eps, dp_air=dp_air, dp_gas=dp_gas, **{k: v for k, v in kw.items() if k not in ("PR", "T41_C")},
                PR=kw["PR"], T41_C=kw["T41_C"])
    r, kw = _scaled_ref(r, kw, P_target)
    k2 = dict(kw); k2.update(p_exh_ratio=1.02 * (1 + dp_gas), dp_rec_air=dp_air, T3r=r["T3r_C"] + 273.15)
    ch = chain(**k2)
    acc = account(ch, env, f3["eta_gb"], f3["eta_gen"], recup=dict(dp_gas=dp_gas, dp_air=dp_air),
                  split=compressor_split(ch, f3["pi1"], f3["eta_p1"]))
    acc["model_check"] = compare_with_model(ch, r)
    acc["inputs"] = dict(kw, eps=eps, dp_air=dp_air, dp_gas=dp_gas)
    acc["cycle"] = dict(ch["out"], T3r_C=r["T3r_C"], Tstack_model_C=r["Tstack_C"])
    return acc, ch


def scan_recup(kw_fn, eps, Tmax, lo=3.5, hi=8.0, step=0.5):
    """Перебор π цикла с рекуператором: грубо с шагом step, затем ±step/2 около лучшей; лучшая — наибольший КПД
    при газе на входе в рекуператор не выше Tmax. kw_fn(π) -> аргументы цикла (охлаждение может зависеть от π)."""
    rows = {}

    def one(pr):
        pr = round(pr, 4)
        if pr in rows:
            return
        kw = kw_fn(pr)
        try:
            r = C.recup(pr, kw["T41_C"], eps, **{k: v for k, v in kw.items() if k not in ("PR", "T41_C")})
        except ValueError:
            rows[pr] = dict(PR=pr, ok=False); return
        rows[pr] = dict(PR=pr, eff=r["eff"], EGT_C=r["EGT_C"], T3r_C=r["T3r_C"], sp_work_kJkg=r["sp_work_kJkg"],
                        x_ngv=kw["x_ngv"], x_rot=kw["x_rot"], ok=r["EGT_C"] <= Tmax)
    n = int(round((hi - lo) / step))
    for i in range(n + 1):
        one(lo + i * step)
    ok = [x for x in rows.values() if x.get("ok")]
    b = max(ok, key=lambda x: x["eff"])["PR"]
    one(b - step / 2); one(b + step / 2)
    out = [rows[k] for k in sorted(rows)]
    ok = [x for x in out if x.get("ok")]
    best = max(ok, key=lambda x: x["eff"]) if ok else None
    return out, best


def scale_W1(acc, P_target_kW):
    return acc["inputs"]["W1"] * P_target_kW / acc["P_shaft_kW"]


def pnp(acc, f3):
    """Проверка «Plug&Play» пакета (мягкий допуск frame3): расход на требуемую мощность и выхлоп."""
    W_req = scale_W1(acc, f3["P_req_kW"])
    T_ex = acc["recuperator"]["T_gas_out_C"] if "recuperator" in acc else acc["cycle"]["EGT_C"]
    return dict(W1_for_P_req=W_req, G_cap=f3["G_cap"], flow_ok=W_req <= f3["G_cap"], T_exhaust_to_package_C=T_ex,
                EGT_cap=f3["EGT_cap"], egt_ok=T_ex <= f3["EGT_cap"],
                P_shaft_at_G_cap_kW=acc["P_shaft_kW"] * f3["G_cap"] / acc["inputs"]["W1"])


def eq_eta_pc(pi1, e1, pi2, e2):
    """Эквивалентный политропный КПД сжатия: та же работа двух ступеней (Cantera)."""
    a = C.gas_air()
    w = C.poly_process(a, P0 * pi1, e1, True) + C.poly_process(a, P0 * pi1 * pi2, e2, True)

    def r(eta):
        b = C.gas_air()
        return C.poly_process(b, P0 * pi1 * pi2, eta, True) - w
    return brentq(r, 0.6, 0.99, xtol=1e-9)


# ---------------------------------------------------------------- SVG: диаграмма Грассмана и сравнение узлов
COL = dict(flow="#3987e5", fuel="#0d366b", work="#6250d6", destr="#eb6834", out="#7d7c78", block="#ecebe7",
           ink="#0b0b0b", ink2="#52514e", grid="#d6d5d0")


HALO = 'stroke="#ffffff" stroke-width="3" paint-order="stroke"'


def _f(x):
    return f"{x:.1f}"


def _pct(x):
    return f"{x:.1f}".replace(".", ",") + " %"


def grassmann_svg(acc, title, subtitle, scale, path):
    """Диаграмма Грассмана: ширина полосы ∝ потоку эксергии (кВт). Слева направо — тракт: компрессор,
    [рекуператор, холодная сторона], камера, турбина, [рекуператор, горячая сторона], выхлоп. Сверху — работа
    (петля на компрессор и полезная мощность), снизу — разрушение эксергии по узлам."""
    n = acc["nodes"]; E = acc["Ex_kW"]; Wk = acc["works_kW"]; Ef = acc["fuel"]["Ex_kW"]
    s = scale
    has_rec = "recuperator" in acc
    has_cool = E["cool"] > 1e-6
    kw = lambda key: n[key]["kW"] if key in n else 0.0  # noqa: E731
    W, Hh = 1120, 660
    yr = 360.0                                   # ось тракта
    gap = 4.0
    el = []
    labels = []                                  # (x_band, y_band_end, text1, text2, kind)

    def band_v(x, y0, y1, w, color, tip, op=0.9):
        if w <= 0:
            return
        el.append(f'<path d="M{_f(x)},{_f(y0)} V{_f(y1)}" stroke="{color}" stroke-width="{max(w, 0.8):.2f}" fill="none" '
                  f'stroke-opacity="{op}"><title>{tip}</title></path>')

    def band_poly(pts, w, color, tip, op=0.9):
        if w <= 0:
            return
        d = "M" + " L".join(f"{_f(x)},{_f(y)}" for x, y in pts)
        el.append(f'<path d="{d}" stroke="{color}" stroke-width="{max(w, 0.8):.2f}" fill="none" stroke-linejoin="miter" '
                  f'stroke-opacity="{op}"><title>{tip}</title></path>')

    def river(x0, x1, flow, tip):
        h = s * flow
        el.append(f'<rect x="{_f(x0)}" y="{_f(yr - h / 2)}" width="{_f(x1 - x0)}" height="{_f(max(h, 0.8))}" fill="{COL["flow"]}" '
                  f'fill-opacity="0.9"><title>{tip}: {flow:.0f} кВт</title></rect>')

    def block(x0, wdt, flow_max, name):
        h = s * flow_max + 16
        el.append(f'<rect x="{_f(x0)}" y="{_f(yr - h / 2)}" width="{_f(wdt)}" height="{_f(h)}" rx="4" fill="{COL["block"]}" '
                  f'stroke="{COL["ink2"]}" stroke-width="1"/>')
        el.append(f'<text x="{_f(x0 + wdt / 2)}" y="{_f(yr + 4)}" text-anchor="middle" font-size="12" fill="{COL["ink"]}">{name}</text>')
        return yr - h / 2, yr + h / 2

    def drops(x0, wdt, ybot, items):
        """полосы разрушения вниз из блока: items = [(kW, подпись, вид)]"""
        widths = [s * v for v, _, _ in items]
        tot = sum(max(w_, 0.8) for w_ in widths) + gap * (len(items) - 1)
        x = x0 + (wdt - tot) / 2
        for (v, lab, kind), w_ in zip(items, widths):
            w2 = max(w_, 0.8)
            xc = x + w2 / 2
            col = COL["destr"] if kind == "d" else COL["out"]
            band_v(xc, ybot, ybot + 46, w_, col, f"{lab}: {v:.0f} кВт, {_pct(100 * v / Ef)} эксергии топлива")
            el.append(f'<path d="M{_f(xc - max(w2 / 2, 3))},{_f(ybot + 46)} L{_f(xc + max(w2 / 2, 3))},{_f(ybot + 46)} '
                      f'L{_f(xc)},{_f(ybot + 54)} Z" fill="{col}"/>')
            labels.append((xc, ybot + 56, lab, _pct(100 * v / Ef)))
            x += w2 + gap

    # --- сначала ширины блоков
    Wc, Wt = Wk["W_comp"], Wk["W_turb"]
    cool = E["cool"]
    air_c = E["air_c3"]                                         # воздух в камеру (или в рекуператор)
    air_in = E["air_comb_in"]
    comb_in = air_in + Ef
    s4, s5, exh = E["s4"], E["s5"], E["exhaust"]
    rec_tr = E["recup_transfer"]
    x = 70.0
    wc_blk = max(s * (Wc + cool) + 30, 90)
    xc0 = x; x += wc_blk + 40
    if has_rec:
        wrc = max(s * rec_tr + 30, 60); xrc0 = x; x += wrc + 40
    comb_items = [(kw("comb_mixing"), "смешение", "d"), (kw("comb_reaction"), "горение", "d"),
                  (kw("comb_dp"), "Δp камеры", "d"), (kw("comb_heat_loss"), "недожог", "d")]
    wb = max(s * Ef + 30, sum(max(s * v, 0.8) for v, _, _ in comb_items) + gap * len(comb_items) + 30, 170)
    xb0 = x; x += wb + 50
    turb_items = [(kw("turbine_gg") + kw("turbine_pt"), "турбина", "d")]
    if has_cool:
        turb_items.append((kw("cool_throttle_ngv") + kw("cool_mix_ngv"), "охл. СА", "d"))
        turb_items.append((kw("cool_throttle_rotor") + kw("cool_mix_rotor"), "охл. ротора", "d"))
    if kw("duct_dp") > 0.5:
        turb_items.append((kw("duct_dp"), "канал", "d"))
    wt = max(s * (Wt + cool) + 30, 120)
    xt0 = x; x += wt + 40
    if has_rec:
        wrg = max(s * rec_tr + 40, 90); xrg0 = x; x += wrg + 40
    x_end = W - 40
    # --- вертикальная раскладка сверху вниз: петля работы, охлаждающий воздух, подпись топлива, тракт
    net = Wt - Wc
    y_loop = 80 + s * Wc / 2
    loop_bot = y_loop + s * Wc / 2
    y_cool = loop_bot + 24 + s * cool / 2
    fuel_top = (y_cool + s * cool / 2 if has_cool else loop_bot) + 46
    yr = fuel_top + 36 + (s * max(comb_in, s4) + 16) / 2
    # --- тракт
    el.append(f'<text x="{_f(xc0 - 52)}" y="{_f(yr - 6)}" font-size="11" fill="{COL["ink2"]}">воздух</text>')
    el.append(f'<text x="{_f(xc0 - 52)}" y="{_f(yr + 8)}" font-size="11" fill="{COL["ink2"]}">{_pct(acc["totals"]["air_in_pct"])}</text>')
    river(xc0 + wc_blk, (xrc0 if has_rec else xb0), air_c, "сжатый воздух в камеру")
    if has_rec:
        river(xrc0 + wrc, xb0, air_in, "подогретый воздух в камеру")
    river(xb0 + wb, xt0, s4, "газ за камерой")
    river(xt0 + wt, (xrg0 if has_rec else x_end - 14), s5, "газ за турбиной")
    if has_rec:
        river(xrg0 + wrg, x_end - 14, exh, "газ за рекуператором")
    hx = s * exh
    el.append(f'<path d="M{_f(x_end - 14)},{_f(yr - hx / 2 - 4)} L{_f(x_end + 4)},{_f(yr)} L{_f(x_end - 14)},{_f(yr + hx / 2 + 4)} Z" fill="{COL["out"]}"/>')
    el.append(f'<text x="{_f(x_end + 6)}" y="{_f(yr - hx / 2 - 22)}" font-size="11" text-anchor="end" fill="{COL["ink"]}">выхлоп</text>')
    el.append(f'<text x="{_f(x_end + 6)}" y="{_f(yr - hx / 2 - 8)}" font-size="11" text-anchor="end" fill="{COL["ink"]}">{_pct(n["exhaust"]["pct_fuel_ex"])}</text>')
    # --- блоки
    yt_c, yb_c = block(xc0, wc_blk, max(air_c + cool + kw("leakage") + kw("compressor"), 1), "компрессор")
    if has_rec:
        yt_rc, yb_rc = block(xrc0, wrc, air_in, "рекуп.")
    yt_b, yb_b = block(xb0, wb, max(comb_in, s4), "камера")
    yt_t, yb_t = block(xt0, wt, max(s4 + cool, s5 + Wt), "турбина")
    if has_rec:
        yt_rg, yb_rg = block(xrg0, wrg, s5, "рекуп.")
    # --- топливо сверху в камеру
    xf = xb0 + wb / 2
    band_v(xf, fuel_top, yt_b, s * Ef, COL["fuel"], f"эксергия топлива: {Ef:.0f} кВт")
    el.append(f'<text x="{_f(xf)}" y="{_f(fuel_top - 22)}" text-anchor="middle" font-size="12" fill="{COL["ink"]}">топливо 100 % = {Ef:.0f} кВт</text>')
    el.append(f'<text x="{_f(xf)}" y="{_f(fuel_top - 8)}" text-anchor="middle" font-size="11" fill="{COL["ink2"]}">(химическая + давление подачи)</text>')
    # --- работа: из турбины вверх; петля на компрессор (влево) и полезная (вправо)
    xw_l = xt0 + 15 + s * cool + s * Wc / 2
    xw_r = xt0 + 15 + s * cool + s * Wc + s * (Wt - Wc) / 2
    xc_in = xc0 + 15 + s * Wc / 2
    band_poly([(xw_l, yt_t), (xw_l, y_loop), (xc_in, y_loop), (xc_in, yt_c)], s * Wc, COL["work"],
              f"работа на привод компрессора: {Wc:.0f} кВт")
    el.append(f'<text x="{_f((xw_l + xc_in) / 2)}" y="{_f(y_loop - s * Wc / 2 - 6)}" text-anchor="middle" font-size="11" '
              f'fill="{COL["ink2"]}">работа на компрессор {Wc:.0f} кВт ({_pct(100 * Wc / Ef)})</text>')
    y_net = max(y_loop + 40, 122 + s * net / 2)
    band_poly([(xw_r, yt_t), (xw_r, y_net), (x_end - 10, y_net)], s * net, COL["work"], f"работа турбины сверх компрессора: {net:.0f} кВт")
    # потери на валу вверх с полосы полезной работы
    xs = xw_r + s * net / 2 + 30
    for key, lab in (("mech", "механика"), ("gearbox", "редуктор"), ("generator", "генератор")):
        v = kw(key); w_ = s * v
        yb = y_net - s * net / 2
        band_v(xs, yb, yb - 34, w_, COL["destr"], f"{lab}: {v:.0f} кВт, {_pct(100 * v / Ef)}")
        el.append(f'<text x="{_f(xs)}" y="{_f(yb - 54)}" text-anchor="middle" font-size="10.5" fill="{COL["ink"]}">{lab}</text>')
        el.append(f'<text x="{_f(xs)}" y="{_f(yb - 41)}" text-anchor="middle" font-size="10.5" fill="{COL["ink2"]}">{_pct(100 * v / Ef)}</text>')
        xs += 70
    el.append(f'<path d="M{_f(x_end - 10)},{_f(y_net - s * net / 2 - 4)} L{_f(x_end + 6)},{_f(y_net)} L{_f(x_end - 10)},{_f(y_net + s * net / 2 + 4)} Z" fill="{COL["work"]}"/>')
    el.append(f'<text x="{_f(x_end + 6)}" y="{_f(y_net + s * net / 2 + 18)}" font-size="12" text-anchor="end" fill="{COL["ink"]}">'
              f'на клеммах {acc["P_el_kW"]:.0f} кВт = {_pct(100 * acc["ex_eff_el"])}</text>')
    # --- охлаждающий воздух: из компрессора поверх тракта в турбину
    if has_cool:
        x_from = xc0 + 15 + s * Wc + 6 + s * cool / 2
        x_to = xt0 + 15 + s * cool / 2
        band_poly([(x_from, yt_c), (x_from, y_cool), (x_to, y_cool), (x_to, yt_t)], s * cool, COL["flow"],
                  f"охлаждающий воздух: {cool:.0f} кВт")
        el.append(f'<text x="{_f((x_from + x_to) / 2)}" y="{_f(y_cool - s * cool / 2 - 5)}" text-anchor="middle" font-size="10.5" '
                  f'fill="{COL["ink2"]}">охлаждающий воздух {_pct(100 * cool / Ef)}</text>')
    # --- разрушение вниз
    drops(xc0, wc_blk, yb_c, [(kw("compressor"), "компрессор", "d"), (kw("leakage"), "утечки", "o")])
    drops(xb0, wb, yb_b, comb_items)
    drops(xt0, wt, yb_t, turb_items)
    if has_rec:
        rec_items = [(kw("recup_heat"), "теплообмен", "d"), (kw("recup_dp_air") + kw("recup_dp_gas"), "Δp рекуп.", "d")]
        drops(xrg0 + 12 + s * rec_tr, wrg - 12 - s * rec_tr, yb_rg, rec_items)
    # подписи разрушения: разнести по горизонтали, выноски
    labels.sort(key=lambda t: t[0])
    xs_l, minsp = [], 66.0
    for i, (xb, yb, a, b) in enumerate(labels):
        xl = xb if not xs_l else max(xb, xs_l[-1] + minsp)
        xs_l.append(xl)
    over = xs_l[-1] - (W - 50) if xs_l else 0
    if over > 0:
        xs_l = [v - over for v in xs_l]
        for i in range(len(xs_l) - 2, -1, -1):
            xs_l[i] = min(xs_l[i], xs_l[i + 1] - minsp)
    ylab = max(t[1] for t in labels) + 26 if labels else yr + 200
    # --- рекуператор: петля снизу (под подписями), от горячей стороны к холодной
    y_bottom = ylab + 13
    if has_rec:
        y_low = ylab + 40 + s * rec_tr / 2
        xg = xrg0 + 12 + s * rec_tr / 2
        xcld = xrc0 + wrc / 2
        band_poly([(xg, yb_rg), (xg, y_low), (xcld, y_low), (xcld, yb_rc)], s * rec_tr, COL["flow"],
                  f"тепло выхлопа воздуху (эксергия): {rec_tr:.0f} кВт")
        el.append(f'<text x="{_f((xg + xcld) / 2)}" y="{_f(y_low + 4)}" text-anchor="middle" font-size="11" '
                  f'fill="{COL["ink"]}">рекуператор возвращает воздуху {rec_tr:.0f} кВт ({_pct(100 * rec_tr / Ef)})</text>')
        y_bottom = y_low + s * rec_tr / 2
    for (xb, yb, a, b), xl in zip(labels, xs_l):
        el.append(f'<path d="M{_f(xb)},{_f(yb)} L{_f(xl)},{_f(ylab - 12)}" stroke="{COL["grid"]}" stroke-width="1" fill="none"/>')
        el.append(f'<text x="{_f(xl)}" y="{_f(ylab)}" text-anchor="middle" font-size="10.5" fill="{COL["ink"]}" {HALO}>{a}</text>')
        el.append(f'<text x="{_f(xl)}" y="{_f(ylab + 13)}" text-anchor="middle" font-size="10.5" fill="{COL["ink2"]}" {HALO}>{b}</text>')
    Hh = int(y_bottom + 50)
    # легенда
    ly = Hh - 18
    lx = 30
    for col, txt in ((COL["fuel"], "топливо"), (COL["flow"], "эксергия потока воздуха и газа"), (COL["work"], "работа"),
                     (COL["destr"], "разрушение эксергии"), (COL["out"], "уходит с потоком")):
        el.append(f'<rect x="{lx}" y="{ly - 9}" width="14" height="10" fill="{col}"/>')
        el.append(f'<text x="{lx + 19}" y="{ly}" font-size="11" fill="{COL["ink2"]}">{txt}</text>')
        lx += 19 + 5.8 * len(txt) + 26
    g = {k: n[k]["pct_fuel_ex"] if k in n else 0.0 for k in n}
    desc = (f"Диаграмма Грассмана, ширина полосы пропорциональна потоку эксергии. Топливо 100 % ({Ef:.0f} кВт). "
            f"Разрушение: горение {_pct(g['comb_reaction'])}, смешение {_pct(g['comb_mixing'])}, компрессор {_pct(g['compressor'])}, "
            f"турбина {_pct(g['turbine_gg'] + g['turbine_pt'])}, механика, редуктор и генератор "
            f"{_pct(g['mech'] + g['gearbox'] + g['generator'])}" + (f", рекуператор {_pct(g['recup_heat'] + g['recup_dp_air'] + g['recup_dp_gas'])}" if has_rec else "")
            + f". Уходит с выхлопом {_pct(g['exhaust'])}, с утечками {_pct(g['leakage'])}. На клеммах {_pct(100 * acc['ex_eff_el'])}.")
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {Hh}" width="{W}" height="{Hh}" role="img" '
           f'aria-labelledby="gt gd" font-family="Times New Roman, Times, serif">\n<title id="gt">{title}</title>\n'
           f'<desc id="gd">{desc}</desc>\n<rect width="{W}" height="{Hh}" fill="#ffffff"/>\n'
           f'<text x="30" y="28" font-size="16" fill="{COL["ink"]}">{title}</text>\n'
           f'<text x="30" y="47" font-size="12" fill="{COL["ink2"]}">{subtitle}</text>\n' + "\n".join(el) + "\n</svg>\n")
    open(path, "w", encoding="utf-8").write(svg)


def bars_svg(variants, path):
    """Сравнение исполнений по группам узлов: доля эксергии топлива, %; три ряда на группу."""
    cols = ["#2a78d6", "#eb6834", "#199e70"]
    keys = [g["key"] for g in variants[0][1]["groups"]] + ["product"]
    labs = [g["label"] for g in variants[0][1]["groups"]] + ["на клеммах (полезно)"]
    vals = []
    for name, acc in variants:
        d = {g["key"]: g["pct_fuel_ex"] for g in acc["groups"]}
        d["product"] = 100 * acc["ex_eff_el"]
        vals.append(d)
    vmax = max(max(v.values()) for v in vals)
    W, row, bh, x0, x1 = 1000, 46, 11, 270, 900
    Hh = 90 + row * len(keys) + 40
    el = []
    sc = (x1 - x0) / (math.ceil(vmax / 5) * 5)
    for t in range(0, int(math.ceil(vmax / 5) * 5) + 1, 5):
        xx = x0 + t * sc
        el.append(f'<line x1="{_f(xx)}" y1="70" x2="{_f(xx)}" y2="{Hh - 40}" stroke="#e4e3df" stroke-width="1"/>')
        el.append(f'<text x="{_f(xx)}" y="{Hh - 24}" text-anchor="middle" font-size="11" fill="#52514e">{t}</text>')
    el.append(f'<text x="{_f((x0 + x1) / 2)}" y="{Hh - 8}" text-anchor="middle" font-size="11" fill="#52514e">% эксергии топлива</text>')
    for i, (k, lab) in enumerate(zip(keys, labs)):
        y = 78 + i * row
        el.append(f'<text x="{x0 - 10}" y="{_f(y + 1.5 * bh + 2)}" text-anchor="end" font-size="12" fill="#0b0b0b">{lab}</text>')
        for j, d in enumerate(vals):
            v = d.get(k, 0.0)
            yy = y + j * (bh + 2)
            wv = max(v * sc, 0.0)
            el.append(f'<rect x="{x0}" y="{_f(yy)}" width="{_f(wv)}" height="{bh}" rx="2" fill="{cols[j]}">'
                      f'<title>{variants[j][0]}: {lab} — {_pct(v)}</title></rect>')
            el.append(f'<text x="{_f(x0 + wv + 4)}" y="{_f(yy + bh - 1.5)}" font-size="10" fill="#52514e">{_pct(v)}</text>')
    lx = 30
    for j, (name, _) in enumerate(variants):
        el.append(f'<rect x="{lx}" y="46" width="14" height="10" fill="{cols[j]}"/>')
        el.append(f'<text x="{lx + 19}" y="55" font-size="12" fill="#0b0b0b">{name}</text>')
        lx += 19 + 6.6 * len(name) + 30
    ce = [d["comb_reaction"] + d["exhaust"] for d in vals]
    title = (f"KV-3: горение и выхлоп — {_pct(min(ce))[:-2]}–{_pct(max(ce))} эксергии топлива во всех трёх исполнениях")
    desc = "Группы узлов, доля эксергии топлива: " + "; ".join(
        f"{lab}: " + ", ".join(f"{variants[j][0]} {_pct(vals[j].get(k, 0.0))}" for j in range(len(vals))) for k, lab in zip(keys, labs))
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {Hh}" width="{W}" height="{Hh}" role="img" '
           f'aria-labelledby="bt bd" font-family="Times New Roman, Times, serif">\n<title id="bt">{title}</title>\n'
           f'<desc id="bd">{desc}</desc>\n<rect width="{W}" height="{Hh}" fill="#ffffff"/>\n'
           f'<text x="30" y="28" font-size="16" fill="#0b0b0b">{title}</text>\n'
           + "\n".join(el) + "\n</svg>\n")
    open(path, "w", encoding="utf-8").write(svg)


# ---------------------------------------------------------------- сборка
def r6(x):
    if isinstance(x, float):
        return round(x, 6)
    if isinstance(x, dict):
        return {k: r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [r6(v) for v in x]
    if isinstance(x, (np.floating,)):
        return round(float(x), 6)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def delta_nodes(a, b, thr=0.05):
    """Узлы, где разрушение выросло/упало больше порога (п. п. эксергии топлива): «сумма прячет потерю»."""
    up, down = [], []
    for k in set(a["nodes"]) | set(b["nodes"]):
        va = a["nodes"].get(k, {}).get("pct_fuel_ex", 0.0)
        vb = b["nodes"].get(k, {}).get("pct_fuel_ex", 0.0)
        if vb - va > thr:
            up.append(dict(node=k, d_pp=vb - va))
        elif va - vb > thr:
            down.append(dict(node=k, d_pp=vb - va))
    return sorted(up, key=lambda x: -x["d_pp"]), sorted(down, key=lambda x: x["d_pp"])


def _num(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",")


def make_svgs(res):
    """Диаграммы из exergy_results.json (тот же словарь, что записан в файл): `python3 exergy.py --svg` перерисовывает
    их без пересчёта."""
    V = res["variants"]
    base, recB, hotC = V["base"], V["recup"], V["hot_printed"]
    eps, T_hot = res["recup_scan"]["eps"], res["hot_sweep"]["chosen_T41_C"]
    variants = [("KV-3 расчётная точка", base), (f"рекуператор ε {_num(eps)}, π {_num(recB['inputs']['PR'], 1)}", recB),
                (f"T41 {T_hot:.0f} °C, печатное охлаждение", hotC)]
    smax = max(max(a["Ex_kW"]["air_comb_in"] + a["Ex_kW"]["fuel"], a["Ex_kW"]["s41"]) for _, a in variants)
    scale = 175.0 / smax
    sub = lambda a: (f"π {_num(a['inputs']['PR'], 1)}, T41 {a['inputs']['T41_C']:.0f} °C, расход {_num(a['inputs']['W1'])} кг/с; "  # noqa: E731
                     f"КПД на клеммах {_pct(100 * a['eff_el'])} (по LHV), эксергетический {_pct(100 * a['ex_eff_el'])}; "
                     f"доли — от эксергии топлива; ISO 15 °C, 60 %")
    tl = lambda a: (f"горение {_pct(a['nodes']['comb_reaction']['pct_fuel_ex'])}, выхлоп {_pct(a['nodes']['exhaust']['pct_fuel_ex'])}, "  # noqa: E731
                    f"на клеммах {_pct(100 * a['ex_eff_el'])} эксергии топлива")
    svgs = res["svg"]
    grassmann_svg(base, "KV-3 (Taurus 60 GS), расчётная точка: " + tl(base), sub(base), scale, os.path.join(HERE, svgs["base"]))
    grassmann_svg(recB, f"KV-3 с рекуператором (ε {_num(eps)}): " + tl(recB), sub(recB), scale, os.path.join(HERE, svgs["recup"]))
    grassmann_svg(hotC, f"KV-3, T41 {T_hot:.0f} °C, печатное охлаждение: " + tl(hotC), sub(hotC), scale,
                  os.path.join(HERE, svgs["hot"]))
    bars_svg(variants, os.path.join(HERE, svgs["bars"]))


def main():
    t_start = time.time()
    f3 = load_frame3()
    env = Env()
    fx = fuel_exergy_check(env)
    print(f"среда ISO: x_H2O={env.x_w:.5f}; эксергия CH4 {fx['ex_ch_kJ_mol_mu']:.2f} кДж/моль "
          f"({fx['ex_ch_MJ_kg']:.3f} МДж/кг, /LHV модели {fx['ratio_ex_LHV_model']:.4f}); "
          f"среда Шаргута {fx['szargut_env_298K_kJ_mol']:.2f}", flush=True)

    # 1. KV-3 в расчётной точке
    kw0 = base_kw(f3)
    base, ch0 = run_simple(kw0, env, f3)
    base["turbine_stages"] = turbine_stages(ch0, f3["turb_ER"], f3["turb_eta"])
    assert abs(base["P_shaft_kW"] - f3["P_shaft_kW_ref"]) < 0.05, (base["P_shaft_kW"], f3["P_shaft_kW_ref"])
    P_eq = base["P_shaft_kW"]
    print(f"KV-3: {base['P_shaft_kW']:.1f} кВт, КПД вал {100 * base['eff_shaft']:.2f} %, эксергетический (клеммы) "
          f"{100 * base['ex_eff_el']:.2f} %, невязка {base['balance']['residual_B_pct_fuel_ex']:.4f} %", flush=True)

    # 2. с рекуператором: та же T41 и КПД узлов, перебор π по наибольшему КПД при газе в рекуператор <= LIM
    eps = 0.85
    Tmax = E.LIM["T_rec_max"]
    rows85, best85 = scan_recup(lambda pr: dict(kw0, PR=pr), eps, Tmax)
    kwB = dict(kw0, PR=best85["PR"])
    recB, chB = run_recup(kwB, env, f3, eps, P_target=P_eq)
    # одна центробежная ступень (π ступени 1 frame3) и ε 0,90 — для таблицы
    recB1, _ = run_recup(dict(kw0, PR=f3["pi1"]), env, f3, eps, P_target=P_eq)
    # компромисс ядра рекуператора (TPMS-решётка компактнее, но с большей потерей давления): ε против Δp
    trade = [dict(eps=eps, dp_air=0.025, dp_gas=0.035, eff=recB["eff_shaft"], EGT_C=recB["cycle"]["EGT_C"],
                  NTU_counterflow_Cr1=eps / (1 - eps))]
    for e_, k_dp in ((0.90, 1.0), (0.85, 1.5), (0.90, 1.5), (0.92, 1.5)):
        dpa, dpg = 0.025 * k_dp, 0.035 * k_dp
        r = C.recup(kwB["PR"], kwB["T41_C"], e_, dp_air=dpa, dp_gas=dpg,
                    **{k: v for k, v in kw0.items() if k not in ("PR", "T41_C")})
        trade.append(dict(eps=e_, dp_air=dpa, dp_gas=dpg, eff=r["eff"], EGT_C=r["EGT_C"],
                          NTU_counterflow_Cr1=e_ / (1 - e_)))
    rec90 = dict(PR=kwB["PR"], eff_shaft=trade[1]["eff"], eff_el=trade[1]["eff"] * f3["eta_gb"] * f3["eta_gen"])
    print(f"рекуператор: π {best85['PR']}, КПД {100 * recB['eff_shaft']:.2f} %; одна ступень π {f3['pi1']}: "
          f"{100 * recB1['eff_shaft']:.2f} %", flush=True)
    print("компромисс ε — Δp: " + "; ".join(f"ε {t['eps']} Δp×{t['dp_gas'] / 0.035:.1f}: {100 * t['eff']:.2f} %" for t in trade), flush=True)

    # 3. повышенная T41 и печатное охлаждение
    xcool, cool_info = cooling_model(f3, base["cycle"]["T3_C"])
    T_hot = 1150.0
    sweep = []
    for T in range(1000, 1251, 25):
        row = dict(T41_C=T)
        for tech in ("cast", "printed"):
            xn, xr = xcool(T, tech)
            r = C.cycle(**dict(kw0, T41_C=float(T), x_ngv=xn, x_rot=xr))
            row[tech] = dict(x_ngv=xn, x_rot=xr, eff=r["eff"], P_kW=r["P_MW"] * 1e3, EGT_C=r["EGT_C"],
                             W1_for_P_eq=kw0["W1"] * P_eq / (r["P_MW"] * 1e3))
        if T in [max(v["T"]) for v in F.COOL.values()] or T == T_hot:
            ck = [k for k, v in F.COOL.items() if T <= max(v["T"])][0] if T <= 1200 else None
            if ck:
                r = C.cycle(**dict(kw0, T41_C=float(T), x_ngv=F.COOL[ck]["xn"], x_rot=F.COOL[ck]["xr"]))
                row["track_class"] = dict(cls=ck, x_ngv=F.COOL[ck]["xn"], x_rot=F.COOL[ck]["xr"], eff=r["eff"],
                                          P_kW=r["P_MW"] * 1e3, EGT_C=r["EGT_C"])
        sweep.append(row)
    print("перебор T41 готов", flush=True)
    xn, xr = xcool(T_hot, "printed")
    kwC = dict(kw0, T41_C=T_hot, x_ngv=xn, x_rot=xr)
    hotC, chC = run_simple(kwC, env, f3, P_target=P_eq)
    # то же T41 с литым охлаждением: по классу K трека и по той же модели при η_c литья — что даёт печать
    hotK, _ = run_simple(dict(kw0, T41_C=T_hot, x_ngv=F.COOL["K"]["xn"], x_rot=F.COOL["K"]["xr"]), env, f3, P_target=P_eq)
    xnc, xrc = xcool(T_hot, "cast")
    hotCc, _ = run_simple(dict(kw0, T41_C=T_hot, x_ngv=xnc, x_rot=xrc), env, f3, P_target=P_eq)
    print(f"T41 {T_hot}: печать {100 * hotC['eff_shaft']:.2f} %, литьё (модель) {100 * hotCc['eff_shaft']:.2f} %, "
          f"класс K трека {100 * hotK['eff_shaft']:.2f} %", flush=True)

    # 4. рекуператор + печатное охлаждение (перебор π при T41 T_hot; доля охлаждения — при T3 этой π)
    def kw_bc(pr):
        a = C.gas_air(); C.poly_process(a, P0 * pr, kw0["eta_pc"], True)
        xn_, xr_ = xcool(T_hot, "printed", T3=a.T - 273.15)
        return dict(kw0, PR=pr, T41_C=T_hot, x_ngv=xn_, x_rot=xr_)
    rowsBC, bestBC = scan_recup(kw_bc, eps, Tmax)
    kwBC = kw_bc(bestBC["PR"])
    recBC, _ = run_recup(kwBC, env, f3, eps, P_target=P_eq)
    print(f"рекуператор + печать: π {bestBC['PR']}, КПД {100 * recBC['eff_shaft']:.2f} %", flush=True)

    # 5. где резерв КПД: приёмы на расчётной точке (тот же расход)
    eta_nofric = eq_eta_pc(f3["pi1"], f3["eta_p1_nofric"], f3["PR"] / f3["pi1"], f3["eta_p2_nofric"])
    pert = [
        ("comp_diffuser", "compressor", dict(eta_pc=eta_nofric),
         "η_p компрессора без трения в узком безлопаточном диффузоре (frame3: eta_basis.d_eta_without_diffuser_friction)",
         "печатный диффузор и колесо, форма — сопряжённая (adjoint) оптимизация канала; шире канал, меньше трения"),
        ("turb_plus1", "turbine", dict(eta_pt_gg=f3["eta_t"] + 0.01, eta_pt_pt=f3["eta_t"] + 0.01),
         "η_p турбины +0,01 (оценка)", "сопряжённая оптимизация 3D-профиля, печатный блиск с бандажом (нет утечки через зазор)"),
        ("leak_half", "leak_cool", dict(x_leak=f3["x_leak"] / 2),
         "утечки 1,5 % -> 0,75 % (оценка)", "цельный печатный корпус: меньше разъёмов, лабиринт и сотовое уплотнение заодно с корпусом"),
        ("dp_comb", "comb_dp", dict(dp_comb=0.025),
         "потеря давления камеры 4 % -> 2,5 % (оценка)", "печатная жаровая труба: подвод воздуха ветвящейся сетью по закону Мюррея, плавный диффузор"),
        ("exh_diffuser", "exhaust", dict(p_exh_ratio=1.01),
         "давление за турбиной 1,02 -> 1,01 p0 (оценка)", "печатный выхлопной диффузор с оптимизированной формой (adjoint)"),
        ("eta_b", "comb_heat_loss", dict(eta_b=0.999), "полнота сгорания 0,995 -> 0,999 (оценка)", "камера: подготовка смеси, время пребывания"),
        ("eta_m", "mech_gen", dict(eta_m=0.995), "механический КПД 0,99 -> 0,995 (оценка)", "один вал, меньше опор; газовые/лепестковые подшипники — проверка нужна"),
        ("fuel_preheat", "comb_mixing", dict(T_fuel=473.15), "подогрев топлива до 200 °C теплом выхлопа (расширение модели)",
         "малый печатный теплообменник топливо — выхлоп (TPMS-гироид)"),
    ]
    reserve = []
    for key, grp, ch_kw, change, how in pert:
        kw1 = dict(kw0); kw1.update(ch_kw)
        acc1, _ = run_simple(kw1, env, f3, split=False)
        g0 = {g["key"]: g["pct_fuel_ex"] for g in base["groups"]}
        g1 = {g["key"]: g["pct_fuel_ex"] for g in acc1["groups"]}
        up, down = delta_nodes(base, acc1)
        reserve.append(dict(key=key, group=grp, change=change, technique=how, params=ch_kw,
                            group_pct_base=g0[grp], group_pct_new=g1[grp],
                            d_eff_shaft_pp=100 * (acc1["eff_shaft"] - base["eff_shaft"]),
                            d_eff_el_pp=100 * (acc1["eff_el"] - base["eff_el"]),
                            d_P_shaft_kW_same_W1=acc1["P_shaft_kW"] - base["P_shaft_kW"],
                            nodes_worse=up, nodes_better=down, balance_pct=acc1["balance"]["residual_B_pct_fuel_ex"]))
        print(f"  {key}: {reserve[-1]['d_eff_shaft_pp']:+.2f} п. п.", flush=True)
    kwall = dict(kw0); [kwall.update(p[2]) for p in pert if p[0] not in ("fuel_preheat",)]
    accall, _ = run_simple(kwall, env, f3, split=False)
    up, down = delta_nodes(base, accall)
    reserve.append(dict(key="all_engine", group="—", change="все приёмы выше, кроме подогрева топлива (новая деталь)",
                        technique="вместе", params={k: v for p in pert if p[0] != "fuel_preheat" for k, v in p[2].items()},
                        d_eff_shaft_pp=100 * (accall["eff_shaft"] - base["eff_shaft"]),
                        d_eff_el_pp=100 * (accall["eff_el"] - base["eff_el"]),
                        d_P_shaft_kW_same_W1=accall["P_shaft_kW"] - base["P_shaft_kW"], nodes_worse=up, nodes_better=down,
                        balance_pct=accall["balance"]["residual_B_pct_fuel_ex"]))
    for key, acc1, change, how in (
            ("recuperator", recB, f"рекуператор ε {eps}, π {kwB['PR']} (перебор, газ в рекуператор <= {E.LIM['T_rec_max']:.0f} °C)",
             "печатный рекуператор с ядром из TPMS-решётки (гироид) переменной плотности; минимум S_gen: баланс ΔT и Δp (Бежан)"),
            ("recuperator_1stage", recB1, f"рекуператор ε {eps}, π {f3['pi1']} — одна центробежная ступень (π ступени 1 frame3)",
             "то же + одна ступень компрессора вместо двух"),
            ("hot_printed", hotC, f"T41 {T_hot:.0f} °C, печатное охлаждение (x_СА {kwC['x_ngv']:.4f}, x_ротора {kwC['x_rot']:.4f})",
             "печатные сопловой аппарат и лопатки с внутренними каналами: ветвление по закону Мюррея, решётки в кромках"),
            ("hot_track_K", hotK, f"T41 {T_hot:.0f} °C, литое охлаждение по классу K трека (5,5 % + 5,5 %)", "сравнение: литьё"),
            ("recup_hot_printed", recBC, f"рекуператор ε {eps} + T41 {T_hot:.0f} °C печатное охлаждение, π {kwBC['PR']}", "вместе")):
        up, down = delta_nodes(base, acc1)
        reserve.append(dict(key=key, group="исполнение", change=change, technique=how,
                            d_eff_shaft_pp=100 * (acc1["eff_shaft"] - base["eff_shaft"]),
                            d_eff_el_pp=100 * (acc1["eff_el"] - base["eff_el"]), d_P_shaft_kW_same_W1=None,
                            nodes_worse=up, nodes_better=down, balance_pct=acc1["balance"]["residual_B_pct_fuel_ex"]))
    reserve.append(dict(key="gearbox_generator", group="mech_gen", change="редуктор и генератор пакета не меняются (Plug&Play)",
                        technique="резерв не снимается двигателем", group_pct_base=base["nodes"]["gearbox"]["pct_fuel_ex"] + base["nodes"]["generator"]["pct_fuel_ex"],
                        d_eff_shaft_pp=0.0, d_eff_el_pp=0.0, d_P_shaft_kW_same_W1=0.0, nodes_worse=[], nodes_better=[]))

    for v in (base, recB, recB1, hotC, hotK, hotCc, recBC):
        v["plug_and_play"] = pnp(v, f3)
        b = v["balance"]
        if not b["closed"]:
            raise SystemExit(f"баланс не закрыт: {b['residual_B_pct_fuel_ex']:.4f} %")
        mc = v["model_check"]
        if mc and (mc["P"] > 1e-6 or mc["Wf"] > 1e-6 or mc["EGT_K"] > 1e-3):
            raise SystemExit(f"цепочка не совпала с gtu_cycle: {mc}")

    svgs = dict(base="grassmann_kv3_base.svg", recup="grassmann_kv3_recup.svg", hot="grassmann_kv3_hot_printed.svg",
                bars="exergy_groups_kv3.svg")

    out = dict(
        title="KV-3: эксергетический анализ — разрушение эксергии по узлам, три исполнения, где резерв КПД",
        method=("0D-цикл: цепочка состояний — копия gtu_cycle.cycle_r (Cantera GRI-Mech 3.0, метан, полное сгорание, политропные "
                "процессы 40 шагами), совпадение с gtu_cycle.cycle/recup проверяется; эксергия потока = тепловая + от давления + "
                "химическая (среда ISO 3977-2: 15 °C, 101,325 кПа, RH 60 %; химические потенциалы среды из Cantera, метан — через "
                "реакцию полного сгорания); разрушение — теорема Гюи — Стодолы T0·S_gen и, независимо, баланс эксергии узла; "
                "охлаждение — x = K·ε/(1-ε) (Holland & Thake 1980; Horlock et al. 2001), K литья — по классу K трека, печать — "
                "K·η_c,литьё/η_c,печать"),
        software=f"python {sys.version.split()[0]}, cantera {ct.__version__}, numpy {np.__version__}, scipy {__import__('scipy').__version__}",
        date=DATE,
        status="оценка уровня облика: КПД узлов, потери, модель охлаждения — принятые величины (оценка); газ — идеальный, "
               "конденсация воды в выхлопе при T0 не учитывается; кинетическая энергия не учитывается (полные параметры)",
        inputs=[
            dict(file=_rel(FRAME3), field="split.pi, inputs[slots…мягкий].T41_C/G/n_out, cycle.frame3.eta_pc/eta_t/dp_duct, "
                 "inputs[gtu_family…].x_leak/eta_gen/eta_gearbox, stage1/stage2.PR/eta_p/eta_basis, turbine.options[3].stages, cycle.required",
                 used={k: v for k, v in f3.items()}),
            dict(file="routes/gtu/gtu_cycle.py", field="cycle, cycle_r, recup, poly_process, AIR, T0, P0, LHV",
                 used=dict(AIR=C.AIR, T0=C.T0, P0=C.P0, LHV=C.LHV)),
            dict(file="routes/gtu/gtu_explore.py", field="LIM['T_rec_max']", used=E.LIM["T_rec_max"]),
            dict(file="routes/gtu/gtu_family.py", field="COOL (классы U, N, K)", used=F.COOL),
            dict(file="routes/gtu/gtu_lineup_results.json", field="lineup[name=KV-3]", used="KV-3 = типоразмер G_frame 25,6 кг/с; расчётная точка — Taurus 60 генератор (frame3)"),
        ],
        assumptions=dict(
            environment=dict(T0_K=env.T, p0_Pa=env.P, RH=RH_ISO, p_sat_Pa=env.p_sat, x_H2O=env.x_w, X=env.X,
                             note="сухой воздух модели (gtu_cycle.AIR) + водяной пар при RH 60 %; p_sat — Cantera ct.Water (15 °C)"),
            fuel_supply="метан при T0 и давлении p3 (как в gtu_cycle: топливный компрессор не считается) — эксергия давления подачи входит в эксергию топлива",
            heat_loss_boundary="недожог η_b модель снимает как теплоту; граница — среда T0: вся эксергия этой теплоты разрушена",
            mech_gen=dict(eta_m="gtu_cycle: 0,99 на обеих частях турбины", eta_gearbox=f3["eta_gb"], eta_gen=f3["eta_gen"],
                          note="редуктор и генератор пакета: KV-3 одновальный, вал 14 950 об/мин = вход редуктора пакета"),
            recuperator=dict(eps=eps, dp_air=0.025, dp_gas=0.035, T_gas_in_max_C=E.LIM["T_rec_max"],
                             note="ε и потери давления — по умолчанию gtu_cycle.recup и перебора трека (оценка); выход газа — "
                                  "из баланса энергии (в модели — оценка через 0,93)"),
            cooling=cool_info,
            equal_power="исполнения пересчитаны на ту же мощность на валу, что и расчётная точка (расход пропорционален)",
            comp_split="ступень 1 — π и η_p frame3, ступень 2 — η_p, при котором выход совпадает с эквивалентным КПД цепочки",
        ),
        fuel_exergy_check=fx,
        variants=dict(base=base, recup=recB, recup_1stage=recB1, hot_printed=hotC, hot_track_K=hotK, hot_cast_model=hotCc,
                      recup_hot_printed=recBC),
        recup_scan=dict(eps=eps, rows=rows85, best=best85, eps090_at_best_PR=rec90, eps_vs_dp=trade,
                        eps_vs_dp_note="NTU противотока при равных теплоёмкостях потоков = ε/(1-ε) (оценка размера ядра: "
                                       "площадь ∝ NTU); Δp×1,5 — плата за компактную решётку (оценка)"),
        recup_hot_scan=dict(T41_C=T_hot, rows=rowsBC, best=bestBC),
        hot_sweep=dict(rows=sweep, chosen_T41_C=T_hot,
                       note="T41 выбрана посередине между классами N и K трека; печатные жаропрочные сплавы — оценка"),
        reserve=reserve,
        svg=svgs,
        feed_to_funnels=dict(
            engine_funnel=dict(
                new_axes=[
                    dict(axis="recuperator_eps", values=[0.0, 0.85, 0.90], note="при ε > 0 оптимум π смещается к 4,5–5 (recup_scan)"),
                    dict(axis="T41_with_cooling_tech", values=[1000, 1050, 1100, 1150, 1200], tech=["cast", "printed"],
                         note="доли охлаждения — cooling_model (оценка); печать — η_c 0,75 против 0,5"),
                    dict(axis="x_leak", values=[0.015, 0.0075], note="цельный печатный корпус (оценка)"),
                    dict(axis="dp_comb", values=[0.04, 0.025], note="печатная жаровая труба (оценка)"),
                    dict(axis="eta_pc_diffuser", values=[f3["eta_pc"], eta_nofric], note="frame3 без трения в узком диффузоре"),
                ],
                new_models=["exergy.account(chain(...)): разрушение по узлам и эксергетический КПД как критерий отбора",
                            "exergy.cooling_model: доля охлаждения от T41, T3 и технологии"],
                criterion="эксергетический КПД на клеммах и наибольшая доля разрушения узла — куда тратить расчёты"),
            combustor_funnel=dict(
                new_axes=[dict(axis="air_preheat_T3r", values="T3 или T3r рекуператора", note="горение — крупнейшее разрушение; подогрев воздуха его снижает"),
                          dict(axis="dp_comb", values=[0.04, 0.03, 0.025]),
                          dict(axis="fuel_preheat_C", values=[15, 200])],
                targets=dict(comb_reaction_pct_base=base["nodes"]["comb_reaction"]["pct_fuel_ex"],
                             comb_mixing_pct_base=base["nodes"]["comb_mixing"]["pct_fuel_ex"],
                             comb_dp_pct_base=base["nodes"]["comb_dp"]["pct_fuel_ex"])),
        ),
        runtime_s=None,
    )
    out["runtime_s"] = round(time.time() - t_start, 1)
    res = r6(out)
    json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    make_svgs(res)
    print(f"готово за {out['runtime_s']} с -> {OUT}")


if __name__ == "__main__":
    if "--svg" in sys.argv[1:]:
        make_svgs(json.load(open(OUT, encoding="utf-8")))
    else:
        main()
