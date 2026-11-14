"""Проверка модели цикла вне калибровки по сверенным паспортам установок (routes/gtu/fleet/*.json).

Модель (gtu_cycle.cycle) откалибрована по двум машинам — SGT-400 13 МВт (привод) и Taurus 60 (привод),
gtu_study.py -> gtu_study_results.json. Здесь при калиброванных КПД узлов и долях охлаждения для каждой
записи паспорта:
  1) вход — опубликованные степень сжатия π, расход и температура выхлопа (условия ISO);
  2) температура газа перед рабочим колесом 1-й ступени T41 подбирается так, чтобы температура выхлопа
     модели совпала с паспортной (brentq);
  3) расход воздуха — так, чтобы расход выхлопа модели совпал с паспортным (модель линейна по расходу:
     удельные величины и КПД от W1 не зависят — test_balances.py);
  4) предсказанные мощность и КПД сравниваются с паспортом; генераторное исполнение — на клеммах
     (модель на валу × КПД генератора и редуктора, принятые, validation_fleet.toml);
  5) «ошибка расхода» — насколько расход выхлопа, нужный модели для паспортной мощности, отличается от
     паспортного: dW = P_паспорт / P_модель − 1 (та же ошибка удельной работы, в единицах расхода).
Рекуперативная установка (Mercury 50) считается и простым циклом (чтобы показать, где модель ломается),
и recup с принятой степенью регенерации; температура за рекуператором — точная, из баланса
(balances.stack_temperature).

Конфигурация — validation_fleet.toml. Результат — validation_fleet_results.json и validation_fleet.csv.
Запуск: python3 validate_fleet.py  (около 5–10 мин на 4 ядрах)."""
from __future__ import annotations

import csv
import datetime
import glob
import json
import math
import sys
import time
import tomllib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import balances as B  # noqa: E402
import gtu_cycle as G  # noqa: E402

CFG_PATH = HERE / "validation_fleet.toml"


def load_cfg(path: Path = CFG_PATH) -> dict:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def rel(p: Path) -> str:
    """Путь от корня репозитория (routes/gtu/...)."""
    p = Path(p).resolve()
    root = HERE.parent.parent.parent
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)


# ---------------------------------------------------------------- калибровка и паспорта
def load_calibration(cfg: dict) -> dict:
    st = json.load(open((HERE / cfg["inputs"]["calibration"]).resolve()))
    sets = {}
    for key in ("T60", "SGT400"):
        s = st[key]
        sets[key] = {k: float(s[k]) for k in ("eta_pc", "eta_pt_gg", "eta_pt_pt", "x_ngv", "x_rot", "x_leak", "T41_C")}
        sets[key]["W_cal"] = float(cfg["calibration_flow_kg_s"][key])
        sets[key]["d"] = float(s["d"])
        sets[key]["name"] = s["name"]
    return sets


def num(x):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    try:
        return float(str(x).replace(",", "").split()[0])
    except Exception:
        return None


def shafts_kind(s) -> str:
    s = str(s or "").lower()
    for k in ("single", "two", "three"):
        if s.startswith(k):
            return k
    return "unknown"


def load_records(cfg: dict) -> list[dict]:
    cal = set(cfg["records"]["calibration"])
    related = set(cfg["records"]["related"])
    recup_models = cfg["records"]["recuperated_models"]
    out = []
    for path in sorted(glob.glob(str((HERE / cfg["inputs"]["fleet_glob"]).resolve()))):
        d = json.load(open(path))
        stem = Path(path).stem.replace("gtu_fleet_", "")
        variants = [str(r.get("variant")) for r in d.get("records", [])]
        for r in d.get("records", []):
            P, eff, W, T, PR = (num(r.get(k)) for k in ("power_kW", "eff_pct", "exh_flow_kg_s", "exh_temp_C", "PR"))
            if not all((P, eff, W, T, PR)):
                continue
            model = r.get("model") or d.get("model")
            rid = f"{stem}:{r.get('variant')}"
            if variants.count(str(r.get("variant"))) > 1:      # SGT-A05: три записи «PG» разных моделей
                rid = f"{rid}:{model}"
            var = str(r.get("variant", "")).upper()
            out.append(dict(
                id=rid, file=rel(path), model=model, variant=r.get("variant"),
                maker="Siemens" if stem.startswith("sgt") else "Solar",
                generator=var.startswith(("GS", "PG")), shafts=shafts_kind(r.get("shafts")),
                shafts_text=r.get("shafts"), P_kW=P, eff_pct=eff, HR_kJ_kWh=num(r.get("heat_rate_kJ_kWh")),
                W_exh=W, EGT_C=T, PR=PR,
                cls="calibration" if rid in cal else ("related" if rid in related else "independent"),
                recuperated=any(m in str(model) for m in recup_models),
                aeroderivative="SGT-A05" in str(model)))
    return out


# ---------------------------------------------------------------- наборы КПД узлов
def interp_t(W: float, cal: dict) -> float:
    a, b = math.log(cal["T60"]["W_cal"]), math.log(cal["SGT400"]["W_cal"])
    return min(1.0, max(0.0, (math.log(W) - a) / (b - a)))


def base_params(rec: dict, var: dict, cal: dict) -> dict:
    keys = ("eta_pc", "eta_pt_gg", "eta_pt_pt", "x_ngv", "x_rot", "x_leak", "T41_C")
    if var["kind"] == "fixed":
        s = cal[var["set"]]
        p = {k: s[k] for k in keys}
        p["t_interp"] = None
    else:
        t = interp_t(rec["W_exh"], cal)
        p = {k: cal["T60"][k] + t * (cal["SGT400"][k] - cal["T60"][k]) for k in keys}
        p["t_interp"] = t
    p["T41_ref_C"] = p.pop("T41_C")
    return p


def cycle_kwargs(p: dict, var: dict, T41: float, rec: dict) -> dict:
    k = 1.0
    if var.get("cool_T41"):
        k = max(0.0, 1 + 0.12 * (T41 - p["T41_ref_C"]) / 50)
    kw = dict(eta_pc=p["eta_pc"], eta_pt_gg=p["eta_pt_gg"], eta_pt_pt=p["eta_pt_pt"],
              x_ngv=p["x_ngv"] * k, x_rot=p["x_rot"] * k, x_leak=p["x_leak"])
    if var.get("single_shaft_no_duct") and rec["shafts"] == "single":
        kw["dp_duct"] = 0.0
    return kw


# ---------------------------------------------------------------- подбор T41 под температуру выхлопа
class NoSolution(Exception):
    pass


def fit_T41(f, target: float, Tmin: float, Tmax: float, xtol: float):
    """f(T41) -> (выход модели, температура выхлопа °C) или None, если точка не считается.
    Температура выхлопа растёт с T41. Возвращает (T41, выход, число вызовов)."""
    cache = {}

    def ev(T):
        T = round(T, 6)
        if T not in cache:
            try:
                cache[T] = f(T)
            except Exception:
                cache[T] = None
        return cache[T]

    a, b = 950.0, 1150.0
    ra = ev(a)
    while ra is None and a < Tmax - 100:
        a += 50
        ra = ev(a)
    if ra is None:
        raise NoSolution("модель не считается ни в одной пробной точке T41")
    b = max(b, a + 100)
    rb = ev(b)
    if rb is None:
        b = a + 50
        rb = ev(b)
    slope = (rb[1] - ra[1]) / (b - a) if rb else 0.55
    T_est = min(Tmax, max(Tmin, a + (target - ra[1]) / slope))
    lo, hi, step = T_est - 30, T_est + 30, 30.0
    for _ in range(20):
        lo, hi = max(Tmin, lo), min(Tmax, hi)
        rl, rh = ev(lo), ev(hi)
        if rl is None:
            # нижняя граница не считается (ТВД не вращает компрессор): двигаемся вверх
            ok = [T for T, r in cache.items() if r is not None]
            if not ok:
                raise NoSolution("модель не считается")
            Tok = min(ok)
            if cache[Tok][1] > target:
                # найти нижнюю считаемую границу между lo и Tok
                l, h = lo, Tok
                for _ in range(12):
                    m = 0.5 * (l + h)
                    if ev(m) is None:
                        l = m
                    else:
                        h = m
                raise NoSolution(
                    f"паспортная температура выхлопа {target:.0f} °C ниже наименьшей достижимой моделью: при T41 ≈ "
                    f"{h:.0f} °C (ниже ТВД не вращает компрессор) выхлоп {cache[round(h, 6)][1]:.0f} °C")
            lo = Tok
            continue
        if rh is None:
            hi = 0.5 * (hi + lo)
            continue
        if rl[1] <= target <= rh[1]:
            break
        if rh[1] < target:
            if hi >= Tmax:
                raise NoSolution(f"температура выхлопа {target:.0f} °C не достигается при T41 до {Tmax:.0f} °C")
            lo, hi, step = hi, hi + 2 * step, 2 * step
        else:
            if lo <= Tmin:
                raise NoSolution(f"температура выхлопа {target:.0f} °C ниже достижимой при T41 от {Tmin:.0f} °C")
            hi, lo, step = lo, lo - 2 * step, 2 * step
    else:
        raise NoSolution("не найден интервал с переменой знака")
    T = brentq(lambda t: ev(t)[1] - target, lo, hi, xtol=xtol)
    r = ev(T)
    return T, r[0], len(cache)


# ---------------------------------------------------------------- одна задача: запись × вариант
def eta_conv(rec: dict, cfg: dict) -> float:
    m = cfg["model"]
    return m["eta_generator"] * m["eta_gearbox"] if rec["generator"] else 1.0


def compare(rec: dict, P_shaft_kW_per_kg_exh: float, eff_shaft: float, cfg: dict) -> dict:
    conv = eta_conv(rec, cfg)
    P_shaft = P_shaft_kW_per_kg_exh * rec["W_exh"]
    P_pred = P_shaft * conv
    eff_pred = eff_shaft * conv * 100
    return dict(P_shaft_kW=P_shaft, P_pred_kW=P_pred, eff_pred_pct=eff_pred, conv=conv,
                dP_pct=(P_pred / rec["P_kW"] - 1) * 100,
                deta_pp=eff_pred - rec["eff_pct"],
                dHR_pct=(rec["eff_pct"] / eff_pred - 1) * 100,
                dW_pct=(rec["P_kW"] / P_pred - 1) * 100,
                W_exh_for_passport_power_kg_s=rec["W_exh"] * rec["P_kW"] / P_pred)


def task_cycle(args):
    rec, var, cal, cfg = args
    t0 = time.time()
    p = base_params(rec, var, cal)
    m = cfg["model"]

    def f(T):
        r = G.cycle(rec["PR"], T, W1=1.0, **cycle_kwargs(p, var, T, rec))
        return r, r["EGT_C"]

    out = dict(id=rec["id"], variant=var["key"], params={k: v for k, v in p.items()})
    try:
        T, r, n = fit_T41(f, rec["EGT_C"], m["T41_min_C"], m["T41_max_C"], m["T41_xtol_K"])
    except NoSolution as e:
        out.update(status="нет решения", reason=str(e), seconds=round(time.time() - t0, 1))
        return out
    kw = cycle_kwargs(p, var, T, rec)
    W1 = rec["W_exh"] / r["W_exh"]
    out.update(status="ok", n_eval=n, T41_C=T, T4_C=r["T4_C"], T3_C=r["T3_C"], T45_C=r["T45_C"],
               EGT_model_C=r["EGT_C"], GG_ER=r["GG_ER"], PT_ER=r["PT_ER"], sp_work_kJ_kg_air=r["sp_work_kJkg"],
               eff_shaft_pct=r["eff"] * 100, W1_kg_s=W1, Wf_kg_s=r["Wf"] * W1, x_ngv_used=kw["x_ngv"],
               x_rot_used=kw["x_rot"], dp_duct_used=kw.get("dp_duct", 0.01),
               **compare(rec, r["P_MW"] * 1e3 / r["W_exh"], r["eff"], cfg))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def task_recup(args):
    rec, var, cal, cfg, eps = args
    t0 = time.time()
    p = base_params(rec, var, cal)
    m = cfg["model"]

    def f(T):
        kw = dict(cycle_kwargs(p, var, T, rec), W1=1.0)
        r = G.recup(rec["PR"], T, eps, **kw)
        Ts = B.stack_temperature(r, kw)
        r["Tstack_exact_C"] = Ts
        return r, Ts

    out = dict(id=rec["id"], variant=f"recup_eps{eps:.2f}", eps=eps, base_variant=var["key"],
               params={k: v for k, v in p.items()})
    try:
        T, r, n = fit_T41(f, rec["EGT_C"], m["T41_min_C"], m["T41_max_C"], m["T41_xtol_K"])
    except NoSolution as e:
        out.update(status="нет решения", reason=str(e), seconds=round(time.time() - t0, 1))
        return out
    W1 = rec["W_exh"] / r["W_exh"]
    out.update(status="ok", n_eval=n, T41_C=T, T3_C=r["T3_C"], T3r_C=r["T3r_C"], EGT_turbine_exit_C=r["EGT_C"],
               Tstack_exact_C=r["Tstack_exact_C"], Tstack_recup_estimate_C=r["Tstack_C"],
               eff_shaft_pct=r["eff"] * 100, W1_kg_s=W1, sp_work_kJ_kg_air=r["sp_work_kJkg"],
               **compare(rec, r["P_MW"] * 1e3 / r["W_exh"], r["eff"], cfg))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def _run(task):
    kind, args = task
    return task_cycle(args) if kind == "cycle" else task_recup(args)


# ---------------------------------------------------------------- статистика
METRICS = ("dP_pct", "deta_pp", "dHR_pct", "dW_pct")


def stats(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    s = dict(n=len(ok), n_failed=len(rows) - len(ok), failed=[r["id"] for r in rows if r.get("status") != "ok"])
    for m in METRICS:
        if not ok:
            continue
        v = np.array([r[m] for r in ok])
        i = int(np.argmax(np.abs(v)))
        s[m] = dict(mean=float(v.mean()), mean_abs=float(np.abs(v).mean()), rms=float(np.sqrt((v ** 2).mean())),
                    max_abs=float(abs(v[i])), max_signed=float(v[i]), max_id=ok[i]["id"],
                    min=float(v.min()), max=float(v.max()))
    return s


def trend(rows: list[dict], x: str, y: str, fx=None) -> dict | None:
    ok = [r for r in rows if r.get("status") == "ok" and r.get(x) is not None]
    if len(ok) < 4:
        return None
    X = np.array([fx(r[x]) if fx else r[x] for r in ok])
    Y = np.array([r[y] for r in ok])
    k, b = np.polyfit(X, Y, 1)
    r2 = 1 - ((Y - (k * X + b)) ** 2).sum() / ((Y - Y.mean()) ** 2).sum()
    return dict(x=x if not fx else f"ln({x})", y=y, n=len(ok), slope=float(k), intercept=float(b), r2=float(r2))


def tags(row: dict, rec: dict, cal: dict) -> list[str]:
    t = []
    if rec["recuperated"]:
        t.append("рекуператор: паспортная температура выхлопа — за рекуператором, простой цикл неприменим")
    t.append({"single": "одновальная (модель двухвальная со свободной СТ)", "two": "двухвальная",
              "three": "трёхвальная (модель двухвальная)"}.get(rec["shafts"], "число валов не сверено"))
    if rec["aeroderivative"]:
        t.append("авиапроизводная (паспорт: Centrax 'Aero-derivative single shaft design')")
    Tlo = min(cal["T60"]["T41_C"], cal["SGT400"]["T41_C"])
    Thi = max(cal["T60"]["T41_C"], cal["SGT400"]["T41_C"])
    if row.get("T41_C") is not None:
        if row["T41_C"] > Thi + 25:
            t.append(f"T41 {row['T41_C']:.0f} °C выше калибровочных ({Tlo:.0f}–{Thi:.0f} °C): доли охлаждения экстраполированы")
        elif row["T41_C"] < Tlo - 75:
            t.append(f"T41 {row['T41_C']:.0f} °C ниже калибровочных ({Tlo:.0f}–{Thi:.0f} °C): доли охлаждения калибровки, "
                     "вероятно, больше нужных")
    PRlo, PRhi = 12.2, 16.8   # π калибровки (gtu_study.py)
    if rec["PR"] < PRlo - 0.5 or rec["PR"] > PRhi + 0.5:
        t.append(f"π {rec['PR']} вне диапазона калибровки {PRlo}–{PRhi}")
    if rec["W_exh"] < cal["T60"]["W_cal"] * 0.9 or rec["W_exh"] > cal["SGT400"]["W_cal"] * 1.1:
        t.append(f"расход {rec['W_exh']} кг/с вне диапазона калибровки {cal['T60']['W_cal']}–{cal['SGT400']['W_cal']} кг/с")
    return t


# ---------------------------------------------------------------- главная
def main():
    t_start = time.time()
    cfg = load_cfg()
    cal = load_calibration(cfg)
    recs = load_records(cfg)
    variants = cfg["variant"]
    primary = next(v["key"] for v in variants if v.get("primary"))
    variants = variants + [dict(key="interp_1shaft_noduct", kind="interp", single_shaft_no_duct=True,
                                label="интерполяция; у одновальных без потери переходного канала СТ (dp_duct = 0)")]
    tasks = []
    for v in variants:
        for rec in recs:
            if v.get("single_shaft_no_duct") and rec["shafts"] != "single":
                continue
            tasks.append(("cycle", (rec, v, cal, cfg)))
    rv = next(v for v in variants if v["key"] == cfg["recup"]["variant"])
    for rec in recs:
        if rec["recuperated"]:
            for eps in cfg["recup"]["eps"]:
                tasks.append(("recup", (rec, rv, cal, cfg, float(eps))))
    # дорогие задачи (recup) — первыми
    tasks.sort(key=lambda t: 0 if t[0] == "recup" else 1)
    with ProcessPoolExecutor(max_workers=int(cfg["model"]["workers"])) as ex:
        results = list(ex.map(_run, tasks))

    by_rec = {r["id"]: dict(r, results={}) for r in recs}
    for res in results:
        by_rec[res["id"]]["results"][res["variant"]] = res
    for rec in by_rec.values():
        for key, row in rec["results"].items():
            row["tags"] = tags(row, rec, cal)

    def rows(vkey, cls=None, pred=None):
        out = []
        for rec in by_rec.values():
            row = rec["results"].get(vkey)
            if row is None or (cls and rec["cls"] not in cls) or (pred and not pred(rec)):
                continue
            out.append(dict(row, **{k: rec[k] for k in ("PR", "W_exh", "generator", "shafts", "maker")}))
        return out

    summary = {}
    for v in variants:
        k = v["key"]
        summary[k] = dict(
            label=v["label"],
            independent=stats(rows(k, ("independent",))),
            independent_no_recup=stats(rows(k, ("independent",), lambda r: not r["recuperated"])),
            related=stats(rows(k, ("related",))),
            calibration=stats(rows(k, ("calibration",))),
        )
    P = primary
    noR = lambda r: not r["recuperated"]  # noqa: E731
    groups = {
        "одновальные": stats(rows(P, ("independent",), lambda r: noR(r) and r["shafts"] == "single")),
        "двухвальные": stats(rows(P, ("independent",), lambda r: r["shafts"] == "two")),
        "трёхвальные": stats(rows(P, ("independent",), lambda r: r["shafts"] == "three")),
        "генераторные": stats(rows(P, ("independent",), lambda r: noR(r) and r["generator"])),
        "приводные": stats(rows(P, ("independent",), lambda r: not r["generator"])),
        "Solar": stats(rows(P, ("independent",), lambda r: noR(r) and r["maker"] == "Solar")),
        "Siemens": stats(rows(P, ("independent",), lambda r: r["maker"] == "Siemens")),
        "рекуперативные (простой цикл)": stats(rows(P, ("independent",), lambda r: r["recuperated"])),
    }
    ind = rows(P, ("independent",), noR)
    trends = [t for t in (
        trend(ind, "T41_C", "deta_pp"), trend(ind, "T41_C", "dP_pct"),
        trend(ind, "PR", "deta_pp"), trend(ind, "PR", "dP_pct"),
        trend(ind, "W_exh", "deta_pp", math.log), trend(ind, "W_exh", "dP_pct", math.log)) if t]
    # охлаждение: та же выборка, вариант с ростом охлаждения по T41
    ind_c = rows("interp_coolT41", ("independent",), noR)
    trends_cool = [t for t in (trend(ind_c, "T41_C", "deta_pp"), trend(ind_c, "T41_C", "dP_pct")) if t]
    # одновальные: эффект переходного канала
    one = {r["id"]: r for r in rows(P, None, lambda r: r["shafts"] == "single" and noR(r))}
    one_nd = {r["id"]: r for r in rows("interp_1shaft_noduct", None, lambda r: noR(r))}
    duct = [dict(id=i, deta_pp_interp=one[i]["deta_pp"], deta_pp_noduct=one_nd[i]["deta_pp"],
                 dP_pct_interp=one[i]["dP_pct"], dP_pct_noduct=one_nd[i]["dP_pct"],
                 change_eff_pp=one_nd[i]["deta_pp"] - one[i]["deta_pp"], change_P_pct=one_nd[i]["dP_pct"] - one[i]["dP_pct"])
            for i in one if i in one_nd and one[i]["status"] == "ok" and one_nd[i]["status"] == "ok"]

    br = cfg["breaks"]
    breaks = []
    for rec in by_rec.values():
        row = rec["results"].get(P)
        if row is None:
            continue
        bad = row["status"] != "ok" or abs(row["deta_pp"]) > br["eff_pp"] or abs(row["dP_pct"]) > br["power_pct"]
        if bad:
            breaks.append(dict(id=rec["id"], cls=rec["cls"], status=row["status"], reason=row.get("reason"),
                               dP_pct=row.get("dP_pct"), deta_pp=row.get("deta_pp"), T41_C=row.get("T41_C"),
                               tags=row["tags"]))

    recup_rows = [res for res in results if res["variant"].startswith("recup_")]

    meta = dict(
        method=("Проверка вне калибровки: модель gtu_cycle.cycle (0D, Cantera GRI-Mech 3.0, метан, ISO) с КПД узлов и "
                "долями охлаждения калибровки gtu_study.py; для каждой записи паспорта T41 подобрана под паспортную "
                "температуру выхлопа (brentq), расход воздуха — под паспортный расход выхлопа; предсказанные мощность и "
                "КПД сравнены с паспортом (генератор — на клеммах через принятые КПД генератора и редуктора). "
                "Рекуператор — gtu_cycle.recup, температура за рекуператором из баланса (balances.stack_temperature)."),
        software=B.SOFTWARE, date=B.DATE,
        inputs=sorted({r["file"] for r in recs}) + [rel(HERE / cfg["inputs"]["calibration"]), rel(CFG_PATH),
                                                     "routes/gtu/gtu_cycle.py", "routes/gtu/validation/balances.py",
                                                     "routes/gtu/validation/validate_fleet.py"],
        units=dict(power="кВт", flow="кг/с", temperature="°C", dP_pct="% к паспорту", deta_pp="п. п. к паспорту",
                   dHR_pct="% расхода тепла к паспорту", dW_pct="% расхода выхлопа для паспортной мощности к паспортному"),
        definitions=dict(
            dP_pct="(мощность модели при паспортных π, расходе и температуре выхлопа) / (паспортная мощность) − 1",
            deta_pp="КПД модели − КПД паспорта (генератор — на клеммах, привод — на валу)",
            dHR_pct="расход тепла модели / паспортный − 1 = КПД паспорта / КПД модели − 1",
            dW_pct="расход выхлопа, при котором модель даёт паспортную мощность, / паспортный − 1",
            cls="calibration — запись вошла в калибровку; related — та же модель в другом исполнении; "
                "independent — паспорт в калибровку не входил"),
        primary_variant=P, runtime_s=round(time.time() - t_start, 1),
    )
    calib_out = {k: {kk: vv for kk, vv in v.items()} for k, v in cal.items()}
    res = dict(**meta, config=cfg, calibration_sets=calib_out, summary=summary, groups_primary_independent=groups,
               trends_primary_independent=trends, trends_coolT41_independent=trends_cool,
               single_shaft_duct_effect=duct, breaks_primary=breaks, recuperated=recup_rows,
               records=list(by_rec.values()))
    json.dump(res, open(HERE / "validation_fleet_results.json", "w"), ensure_ascii=False, indent=1, default=float)
    write_csv(by_rec, variants, primary)
    print_summary(res)


CSV_COLS = ["id", "model", "variant_passport", "cls", "maker", "drive", "shafts", "model_variant", "status",
            "PR", "W_exh_pass_kg_s", "EGT_pass_C", "P_pass_kW", "eff_pass_pct", "T41_fit_C", "EGT_model_C",
            "P_pred_kW", "eff_pred_pct", "dP_pct", "deta_pp", "dHR_pct", "dW_pct", "W_exh_for_passport_power_kg_s",
            "W1_kg_s", "eta_pc", "eta_pt_gg", "eta_pt_pt", "x_ngv", "x_rot", "x_leak", "reason"]


def write_csv(by_rec, variants, primary):
    order = [primary] + [v["key"] for v in variants if v["key"] != primary]
    with open(HERE / "validation_fleet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLS)
        for vk in order + ["recup_eps0.85", "recup_eps0.90"]:
            for rec in by_rec.values():
                row = rec["results"].get(vk)
                if row is None:
                    continue
                p = row["params"]

                def f(x, nd=3):
                    return "" if x is None else (round(x, nd) if isinstance(x, float) else x)
                w.writerow([rec["id"], rec["model"], rec["variant"], rec["cls"], rec["maker"],
                            "генератор" if rec["generator"] else "привод", rec["shafts"], vk, row["status"],
                            rec["PR"], rec["W_exh"], rec["EGT_C"], rec["P_kW"], rec["eff_pct"],
                            f(row.get("T41_C"), 1), f(row.get("EGT_model_C", row.get("Tstack_exact_C")), 1),
                            f(row.get("P_pred_kW"), 0), f(row.get("eff_pred_pct"), 2), f(row.get("dP_pct"), 2),
                            f(row.get("deta_pp"), 2), f(row.get("dHR_pct"), 2), f(row.get("dW_pct"), 2),
                            f(row.get("W_exh_for_passport_power_kg_s"), 2), f(row.get("W1_kg_s"), 3),
                            f(p["eta_pc"], 4), f(p["eta_pt_gg"], 4), f(p["eta_pt_pt"], 4),
                            f(row.get("x_ngv_used", p["x_ngv"]), 4), f(row.get("x_rot_used", p["x_rot"]), 4),
                            f(p["x_leak"], 4), row.get("reason", "")])


def print_summary(res):
    P = res["primary_variant"]
    for rec in res["records"]:
        row = rec["results"].get(P)
        if row["status"] == "ok":
            print(f"{rec['id']:22s} {rec['cls'][:5]} π{rec['PR']:5.1f} G{rec['W_exh']:6.2f} T5 {rec['EGT_C']:4.0f} | "
                  f"T41 {row['T41_C']:6.1f}  P {row['P_pred_kW']:7.0f}/{rec['P_kW']:6.0f} ({row['dP_pct']:+5.1f} %)  "
                  f"η {row['eff_pred_pct']:5.2f}/{rec['eff_pct']:5.2f} ({row['deta_pp']:+5.2f} п.п.)  dW {row['dW_pct']:+5.1f} %")
        else:
            print(f"{rec['id']:22s} {rec['cls'][:5]} — {row['status']}: {row['reason']}")
    for k, s in res["summary"].items():
        i = s["independent_no_recup"]
        if i.get("n"):
            print(f"[{k}] независимые без рекуп. n={i['n']}: |dP| ср {i['dP_pct']['mean_abs']:.2f} макс {i['dP_pct']['max_abs']:.2f} %;"
                  f" |dη| ср {i['deta_pp']['mean_abs']:.2f} макс {i['deta_pp']['max_abs']:.2f} п.п.; смещение dP {i['dP_pct']['mean']:+.2f}, dη {i['deta_pp']['mean']:+.2f}")
    for r in res["recuperated"]:
        if r["status"] == "ok":
            print(f"recup ε={r['eps']}: T41 {r['T41_C']:.0f} P {r['P_pred_kW']:.0f} ({r['dP_pct']:+.1f} %) η {r['eff_pred_pct']:.2f} ({r['deta_pp']:+.2f})")
        else:
            print(f"recup ε={r['eps']}: {r['reason']}")
    print("время", res["runtime_s"], "с")


if __name__ == "__main__":
    main()
