"""Балансы массы и энергии модели цикла трека gas-turbine (routes/gtu/gtu_cycle.py).

Проверяются функции cycle, cycle_r и recup. Модель не меняется: функции вызываются как есть,
а балансы сводятся по их выходу и по входным параметрам.

Граница — весь двигатель: вход воздуха (T0, P0) и топлива (метан, T0), выход выхлопа, утечек
(воздух за компрессором, T3), работы на валу и потерь (механические, недожог).

Баланс массы:   W1 + Wf = W_exh + x_leak·W1.
Баланс энергии двумя способами:
  1) по абсолютным энтальпиям Cantera (с теплотами образования):
     W1·h_air(T0) + Wf·h_CH4(T0) + Q_ext = W_exh·h_exh(T5) + x_leak·W1·h_air(T3) + P + L_mech + L_comb;
  2) в форме «теплота топлива = работа + теплота выхлопа + потери», теплота топлива Wf·LHV
     (LHV = 50,03 МДж/кг — константа модели), теплота выхлопа отсчитана от T0:
     Wf·LHV + Q_ext = P + W_exh·[h_exh(T5) − h_exh(T0)] + x_leak·W1·[h_air(T3) − h_air(T0)] + L_mech + L_comb.
Разница способов 1 и 2 обязана равняться (1 − LHV_Cantera(T0)/LHV)·100 % — несогласованность константы
LHV модели со свойствами Cantera; это вторая, независимая проверка.

Состав выхлопа восстанавливается заново по входу (воздух W1·(1 − x_leak) плюс метан Wf, полное
сгорание CH4 + 2 O2 -> CO2 + 2 H2O), а не берётся из модели: модель сжигает топливо в части воздуха и
потом подмешивает охлаждение, итог по составу должен совпасть.

Q_ext — теплота, подведённая к воздуху перед камерой вне границы (cycle_r с заданной T3r: газовая
сторона рекуператора в cycle_r не моделируется). Для recup рекуператор внутри границы: выхлоп уходит
при температуре за рекуператором; точная температура считается из баланса рекуператора и
сравнивается с оценкой recup (Tstack_C = EGT − 0,93·(T3r − T3)).

Запуск: python3 balances.py  ->  balances_results.json, balances.csv."""
from __future__ import annotations

import csv
import datetime
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import cantera as ct  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402

import gtu_cycle as G  # noqa: E402  (импорт делает chdir в routes/gtu — пути ниже абсолютные)

DATE = datetime.date.today().isoformat()
SOFTWARE = (f"python {sys.version.split()[0]}, cantera {ct.__version__}, numpy {np.__version__}, "
            f"scipy {scipy.__version__}")

_GAS = None


def _gas():
    global _GAS
    if _GAS is None:
        _GAS = ct.Solution("gri30.yaml")
    return _GAS


def air_Y():
    g = _gas()
    g.TPX = G.T0, G.P0, G.AIR
    return g.Y.copy()


def exhaust_Y(W1: float, Wf: float, x_leak: float):
    """Состав выхлопа по входу: воздух W1(1 − x_leak) + метан Wf, полное сгорание."""
    g = _gas()
    i = g.species_index
    MW = g.molecular_weights
    m = W1 * (1 - x_leak) * air_Y()
    n_f = Wf / MW[i("CH4")]
    m[i("O2")] -= 2 * n_f * MW[i("O2")]
    m[i("CO2")] += n_f * MW[i("CO2")]
    m[i("H2O")] += 2 * n_f * MW[i("H2O")]
    if m[i("O2")] <= 0:
        raise ValueError("не хватает кислорода для полного сгорания")
    return m / m.sum()


def h_mass(T: float, Y=None, X=None) -> float:
    """Удельная энтальпия идеального газа, Дж/кг (от давления не зависит)."""
    g = _gas()
    if X is not None:
        g.TPX = T, G.P0, X
    else:
        g.TPY = T, G.P0, Y
    return g.enthalpy_mass


def T_from_h(h: float, Y) -> float:
    g = _gas()
    g.HPY = h, G.P0, Y
    return g.T


def lhv_cantera(T: float = G.T0) -> float:
    """Низшая теплота сгорания метана по Cantera при температуре T, Дж/кг (вода — пар)."""
    g = _gas()
    g.TP = T, G.P0
    hm = g.standard_enthalpies_RT * ct.gas_constant * T  # Дж/кмоль
    i = g.species_index
    dh = hm[i("CH4")] + 2 * hm[i("O2")] - hm[i("CO2")] - 2 * hm[i("H2O")]
    return dh / g.molecular_weights[i("CH4")]


def run_point(kind: str, kw: dict) -> dict:
    """Вызвать функцию модели: kind ∈ cycle / cycle_r / recup; kw — её аргументы."""
    f = {"cycle": G.cycle, "cycle_r": G.cycle_r, "recup": G.recup}[kind]
    return f(**kw)


def balance(kind: str, kw: dict, r: dict | None = None) -> dict:
    """Балансы массы и энергии для одной точки. Возвращает члены (МВт, кг/с) и невязки (%)."""
    if r is None:
        r = run_point(kind, kw)
    W1 = kw.get("W1", 1.0)
    x_ngv, x_rot, x_leak = kw["x_ngv"], kw["x_rot"], kw["x_leak"]
    eta_m = kw.get("eta_m", 0.99)
    eta_b = kw.get("eta_b", 0.995)
    Wf, Wexh, P = r["Wf"], r["W_exh"], r["P_MW"] * 1e6
    T3, T5 = r["T3_C"] + 273.15, r["EGT_C"] + 273.15
    wc = r["wc_kJkg"] * 1e3

    # масса
    m_in, m_out = W1 + Wf, Wexh + x_leak * W1
    mass_res = (m_in - m_out) / m_in * 100

    # энтальпии
    Yx = exhaust_Y(W1, Wf, x_leak)
    h5, h5_0 = h_mass(T5, Y=Yx), h_mass(G.T0, Y=Yx)
    ha0, ha3 = h_mass(G.T0, X=G.AIR), h_mass(T3, X=G.AIR)
    hf0 = h_mass(G.T0, X="CH4:1")
    W_comb_air = W1 * (1 - x_ngv - x_rot - x_leak)

    T3r = None
    if kind == "cycle_r" and kw.get("T3r"):
        T3r = kw["T3r"]
    if kind == "recup":
        T3r = r["T3r_C"] + 273.15
    Q_air_heat = W_comb_air * (h_mass(T3r, X=G.AIR) - ha3) if T3r else 0.0
    Q_ext = Q_air_heat if kind == "cycle_r" else 0.0   # в recup — внутренний перенос

    Q_fuel = Wf * G.LHV
    L_mech = wc * W1 * (1 / eta_m - 1) + P * (1 / eta_m - 1)
    L_comb = (1 - eta_b) * Wf * G.LHV
    Q_leak = x_leak * W1 * (ha3 - ha0)

    out = dict(kind=kind, W1_kg_s=W1, Wf_kg_s=Wf, W_exh_kg_s=Wexh, P_MW=P / 1e6, eff_pct=r["eff"] * 100,
               T3_C=r["T3_C"], EGT_C=r["EGT_C"], mass_residual_pct=mass_res)

    if kind == "recup":
        # точная температура за рекуператором: W_exh·[h(T5) − h(Tstack)] = W_comb_air·[h_air(T3r) − h_air(T3)]
        h_stack = h5 - Q_air_heat / Wexh
        T_stack = T_from_h(h_stack, Yx)
        Q_exh = Wexh * (h_stack - h5_0)
        Q_exh_est = Wexh * (h_mass(r["Tstack_C"] + 273.15, Y=Yx) - h5_0)
        res_abs = (W1 * ha0 + Wf * hf0) - (Wexh * h_stack + x_leak * W1 * ha3 + P + L_mech + L_comb)
        res_est = Q_fuel - (P + Q_exh_est + Q_leak + L_mech + L_comb)
        out.update(Q_recup_MW=Q_air_heat / 1e6, T3r_C=r["T3r_C"], Tstack_exact_C=T_stack - 273.15,
                   Tstack_recup_estimate_C=r["Tstack_C"], Tstack_estimate_error_K=r["Tstack_C"] - (T_stack - 273.15),
                   energy_residual_recup_estimate_pct=res_est / Q_fuel * 100)
    else:
        Q_exh = Wexh * (h5 - h5_0)
        res_abs = (W1 * ha0 + Wf * hf0 + Q_ext) - (Wexh * h5 + x_leak * W1 * ha3 + P + L_mech + L_comb)

    res_lhv = Q_fuel + Q_ext - (P + Q_exh + Q_leak + L_mech + L_comb)
    lhv_ct = lhv_cantera(G.T0)
    out.update(Q_fuel_MW=Q_fuel / 1e6, Q_ext_MW=Q_ext / 1e6, Q_exh_MW=Q_exh / 1e6, Q_leak_MW=Q_leak / 1e6,
               L_mech_MW=L_mech / 1e6, L_comb_MW=L_comb / 1e6,
               energy_residual_abs_pct=res_abs / Q_fuel * 100,
               energy_residual_lhv_pct=res_lhv / Q_fuel * 100,
               lhv_model_MJ_kg=G.LHV / 1e6, lhv_cantera_T0_MJ_kg=lhv_ct / 1e6,
               lhv_form_expected_pct=(1 - lhv_ct / G.LHV) * 100)
    return out


def stack_temperature(r: dict, kw: dict) -> float:
    """Точная температура выхлопа за рекуператором, °C, для результата recup (баланс рекуператора)."""
    W1 = kw.get("W1", 1.0)
    Wf, Wexh = r["Wf"], r["W_exh"]
    Yx = exhaust_Y(W1, Wf, kw["x_leak"])
    T3, T3r, T5 = r["T3_C"] + 273.15, r["T3r_C"] + 273.15, r["EGT_C"] + 273.15
    W_comb_air = W1 * (1 - kw["x_ngv"] - kw["x_rot"] - kw["x_leak"])
    Q = W_comb_air * (h_mass(T3r, X=G.AIR) - h_mass(T3, X=G.AIR))
    return T_from_h(h_mass(T5, Y=Yx) - Q / Wexh, Yx) - 273.15


# --- точки проверки: калибровочные (gtu_study_results.json), расчётная точка типоразмера № 3, крайние ---
def points() -> list[dict]:
    st = json.load(open(HERE.parent / "gtu_study_results.json"))
    pts = []
    for key, W_exh_pass, divisor in (("SGT400", 39.4, 1.006), ("T60", 21.63, 1.004)):
        s = st[key]
        pts.append(dict(id=f"cycle:{key}-калибровка", kind="cycle",
                        source=f"gtu_study_results.json -> {key}; W1 = {W_exh_pass}/{divisor} как в gtu_study.py",
                        kw=dict(PR=s["PR"], T41_C=s["T41_C"], eta_pc=s["eta_pc"], eta_pt_gg=s["eta_pt_gg"],
                                eta_pt_pt=s["eta_pt_pt"], x_ngv=s["x_ngv"], x_rot=s["x_rot"], x_leak=s["x_leak"],
                                W1=W_exh_pass / divisor)))
    pts += [
        dict(id="cycle:типоразмер3-без-охлаждения", kind="cycle",
             source="README: расчётная точка Taurus 60 генератор, π 13,5, 1000 °C, 25,6 кг/с; КПД узлов как в gtu_family.py (2В-2Ц)",
             kw=dict(PR=13.5, T41_C=1000.0, eta_pc=0.860, eta_pt_gg=0.865, eta_pt_pt=0.885,
                     x_ngv=0.0, x_rot=0.0, x_leak=0.015, W1=25.6)),
        dict(id="cycle:малая-π-охлаждаемый-СА", kind="cycle",
             source="навык gtu-cycle, шаг 1: π 6,75, 1050 °C, класс N",
             kw=dict(PR=6.75, T41_C=1050.0, eta_pc=0.855, eta_pt_gg=0.865, eta_pt_pt=0.885,
                     x_ngv=0.03, x_rot=0.0, x_leak=0.015, W1=1.0)),
        dict(id="cycle:большие-потери", kind="cycle",
             source="крайняя точка: π 20, 1250 °C, класс K, недожог 2 %, мех. КПД 0,97, потери давления увеличены",
             kw=dict(PR=20.0, T41_C=1250.0, eta_pc=0.88, eta_pt_gg=0.87, eta_pt_pt=0.89,
                     x_ngv=0.055, x_rot=0.055, x_leak=0.03, eta_b=0.98, eta_m=0.97, dp_comb=0.06, dp_duct=0.02,
                     p_exh_ratio=1.04, W1=40.0)),
        dict(id="cycle_r:подогрев-воздуха-200K", kind="cycle_r",
             source="cycle_r с внешним подогревом воздуха перед камерой: T3r = T3 + 200 K, потеря давления 2,5 %",
             kw=dict(PR=10.0, T41_C=1100.0, eta_pc=0.874, eta_pt_gg=0.864, eta_pt_pt=0.884,
                     x_ngv=0.045, x_rot=0.045, x_leak=0.015, W1=17.9, dp_rec_air=0.025, T3r="T3+200")),
        dict(id="recup:Mercury-подобный", kind="recup",
             source="recup: π 9,9, 1150 °C, ε 0,90 (принято), КПД узлов Taurus 60 (gtu_study_results.json -> T60)",
             kw=dict(PR=9.9, T41_C=1150.0, eps=0.90, eta_pc=st["T60"]["eta_pc"], eta_pt_gg=st["T60"]["eta_pt_gg"],
                     eta_pt_pt=st["T60"]["eta_pt_pt"], x_ngv=0.045, x_rot=0.045, x_leak=0.015, W1=17.9)),
        dict(id="recup:схема-130-подобный", kind="recup",
             source="recup: π 4,25, 950 °C, ε 0,85 (как схема 130 в gtu_upstream.py), без охлаждения",
             kw=dict(PR=4.25, T41_C=950.0, eps=0.85, eta_pc=0.855, eta_pt_gg=0.865, eta_pt_pt=0.885,
                     x_ngv=0.0, x_rot=0.0, x_leak=0.015, W1=10.0)),
    ]
    return pts


def resolve_kw(p: dict) -> dict:
    """«T3+200» -> температура за компрессором этой же точки плюс 200 K (считается cycle без подогрева)."""
    kw = dict(p["kw"])
    if isinstance(kw.get("T3r"), str):
        base = {k: v for k, v in kw.items() if k not in ("T3r", "dp_rec_air")}
        T3 = G.cycle(**base)["T3_C"] + 273.15
        kw["T3r"] = T3 + float(kw["T3r"].split("+")[1])
    return kw


def main():
    rows = []
    for p in points():
        kw = resolve_kw(p)
        b = balance(p["kind"], kw)
        b.update(id=p["id"], source=p["source"], inputs={k: v for k, v in kw.items()})
        rows.append(b)
        print(f"{p['id']:34s} масса {b['mass_residual_pct']:+.2e} %  энергия(абс) {b['energy_residual_abs_pct']:+.2e} %"
              f"  энергия(LHV) {b['energy_residual_lhv_pct']:+.4f} %"
              + (f"  Tstack точн. {b['Tstack_exact_C']:.1f} / оценка {b['Tstack_recup_estimate_C']:.1f} °C"
                 f"  невязка с оценкой {b['energy_residual_recup_estimate_pct']:+.3f} %" if p["kind"] == "recup" else ""))
    summary = dict(
        max_abs_mass_residual_pct=max(abs(r["mass_residual_pct"]) for r in rows),
        max_abs_energy_residual_abs_pct=max(abs(r["energy_residual_abs_pct"]) for r in rows),
        max_abs_energy_residual_lhv_pct=max(abs(r["energy_residual_lhv_pct"]) for r in rows),
        max_abs_lhv_form_minus_expected_pct=max(abs(r["energy_residual_lhv_pct"] - r["energy_residual_abs_pct"]
                                                    - r["lhv_form_expected_pct"]) for r in rows),
        recup_Tstack_estimate_error_K=[r["Tstack_estimate_error_K"] for r in rows if r["kind"] == "recup"],
        recup_energy_residual_with_estimate_pct=[r["energy_residual_recup_estimate_pct"] for r in rows if r["kind"] == "recup"],
    )
    meta = dict(
        method=("Балансы массы и энергии по выходу gtu_cycle.cycle / cycle_r / recup (модель без изменений); "
                "энергия двумя способами: абсолютные энтальпии Cantera и форма «теплота топлива = работа + теплота "
                "выхлопа + потери» с LHV модели; состав выхлопа восстановлен по входу (полное сгорание)"),
        software=SOFTWARE, date=DATE,
        inputs=["routes/gtu/gtu_cycle.py", "routes/gtu/gtu_study_results.json", "routes/gtu/validation/balances.py"],
        units=dict(mass="кг/с", energy="МВт", residual="% от расхода массы на входе / от теплоты топлива Wf·LHV",
                   temperature="°C"),
    )
    json.dump(dict(**meta, summary=summary, points=rows), open(HERE / "balances_results.json", "w"),
              ensure_ascii=False, indent=1)
    cols = ["id", "kind", "W1_kg_s", "Wf_kg_s", "W_exh_kg_s", "P_MW", "eff_pct", "Q_fuel_MW", "Q_ext_MW", "Q_exh_MW",
            "Q_leak_MW", "L_mech_MW", "L_comb_MW", "mass_residual_pct", "energy_residual_abs_pct",
            "energy_residual_lhv_pct", "lhv_form_expected_pct", "Tstack_exact_C", "Tstack_recup_estimate_C",
            "Tstack_estimate_error_K", "energy_residual_recup_estimate_pct"]
    with open(HERE / "balances.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([("" if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r.get(c), float) else r[c])) for c in cols])
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
