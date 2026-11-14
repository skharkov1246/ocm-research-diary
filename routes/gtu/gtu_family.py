"""Семейство двигателей-заменителей для линейки Solar и Siemens до 20 МВт.
Критерий владельца (14.11.2026): минимум движущихся частей, безусловная надёжность; мощность, частота
выходного вала и габарит - как у исходной установки; КПД может быть ниже.

Для каждой установки (модель x исполнение) перебираются архитектуры по возрастанию сложности:
  1В-1Ц  одновальный, одно центробежное колесо + турбина (только генератор: частота вала задана редуктором);
  2В-1Ц  двухвальный: газогенератор (одно центробежное колесо + одна ступень турбины) + свободная силовая турбина;
  2В-2Ц  двухвальный, два центробежных колеса на валу газогенератора;
и классы охлаждения: U - без охлаждения (T41 <= 1000 C), N - охлаждаемый только сопловой аппарат
(неподвижная литая деталь, T41 <= 1050 C), K - охлаждаемые сопловой аппарат и рабочая лопатка 1-й ступени (T41 <= 1200 C).
Для допуска по расходу alpha (прирост расхода к исходному) ищется наименьшая температура газа, при которой
удельная работа даёт исходную мощность; при равенстве - наибольший КПД.
Модель цикла - gtu_cycle.cycle (калибровка по SGT-400 и Taurus 60). Одновальная схема считается тем же циклом
(без потерь переходного канала - оценка в пользу двухвальной не более 0,3 п. п.).
Вход: файлы паспортов (аргументы командной строки или fleet/*.json). Выход: gtu_family_results.json."""
import json, math, glob, sys, os
from concurrent.futures import ProcessPoolExecutor
from gtu_cycle import cycle
import gtu_explore as E
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))

ARCH = {
    "1В-1Ц": dict(name="одновальный, одно центробежное колесо", rotors=1, nc=1, eta_pc=0.855, PR=[x / 4 for x in range(12, 28)]),
    "1В-2Ц": dict(name="одновальный, два центробежных колеса", rotors=1, nc=2, eta_pc=0.860, PR=[x / 2 for x in range(12, 33)]),
    "2В-1Ц": dict(name="двухвальный, одно центробежное колесо", rotors=2, nc=1, eta_pc=0.855, PR=[x / 4 for x in range(12, 28)]),
    "2В-2Ц": dict(name="двухвальный, два центробежных колеса", rotors=2, nc=2, eta_pc=0.860, PR=[x / 2 for x in range(12, 33)]),
}
COOL = {
    "U": dict(name="без охлаждения", T=[900, 925, 950, 975, 1000], xn=0.0, xr=0.0, rank=0),
    "N": dict(name="охлаждаемый сопловой аппарат", T=[1025, 1050], xn=0.03, xr=0.0, rank=1),
    "K": dict(name="охлаждаемые сопловой аппарат и рабочая лопатка", T=[1075, 1100, 1125, 1150, 1175, 1200], xn=0.055, xr=0.055, rank=2),
}
SCEN = {   # допуск: прирост расхода через вход и выхлоп пакета; прирост температуры выхлопа
    "строгий": dict(alpha=0.0, dT=25),
    "базовый": dict(alpha=0.10, dT=50),
    "мягкий": dict(alpha=0.20, dT=100),
}
D2_MAX = 0.90                    # м, наибольшая поковка титанового колеса (принято)
RECUP = ("Mercury",)             # рекуперативные установки: рекуператор пакета сохраняется, считаются отдельно
ETA_GEN = 0.972 * 0.985          # генератор и редуктор: перевод кВт на клеммах в кВт на валу двигателя
FAR = 0.017                      # доля топлива в расходе выхлопа (оценка): расход воздуха = выхлоп / (1 + FAR)


def _run(args):
    ak, ck, PR, T = args
    a, c = ARCH[ak], COOL[ck]
    try:
        r = cycle(PR, T, eta_pc=a["eta_pc"], eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=c["xn"], x_rot=c["xr"], x_leak=0.015, W1=1.0)
        return dict(arch=ak, cool=ck, PR=PR, T41=T, w=r["sp_work_kJkg"], eta=r["eff"] * 100, EGT=r["EGT_C"], T3=r["T3_C"],
                    GG_ER=r["GG_ER"], PT_ER=r["PT_ER"])
    except Exception:
        return None


def grid():
    jobs = []
    for ak in ("1В-1Ц", "2В-2Ц"):          # 2В-1Ц имеет тот же цикл, что 1В-1Ц
        for ck, c in COOL.items():
            for PR in ARCH[ak]["PR"]:
                for T in c["T"]:
                    jobs.append((ak, ck, PR, T))
    with ProcessPoolExecutor(max_workers=4) as ex:
        pts = [p for p in ex.map(_run, jobs, chunksize=8) if p]
    pts += [dict(p, arch="2В-1Ц") for p in pts if p["arch"] == "1В-1Ц"]
    pts += [dict(p, arch="1В-2Ц") for p in pts if p["arch"] == "2В-2Ц"]
    return pts


def num(x):
    try:
        return float(str(x).replace(",", "").split()[0])
    except Exception:
        return None


def rpm_of(rec):
    """номинальная частота выходного вала: первое число в поле out_rpm"""
    v = rec.get("out_rpm")
    if isinstance(v, (int, float)):
        return float(v)
    if not v:
        return None
    import re
    m = re.findall(r"\d[\d\s,]*", str(v))
    vals = [float(x.replace(",", "").replace(" ", "")) for x in m if x.replace(",", "").replace(" ", "").isdigit()]
    vals = [x for x in vals if 1000 <= x <= 40000]
    return vals[0] if vals else None


def stages_sum(v):
    import re
    if v is None:
        return 0
    if isinstance(v, (int, float)):
        return int(v)
    return sum(int(x) for x in re.findall(r"\d+", str(v)) if int(x) < 30)


def slots(paths):
    out = []
    for p in paths:
        d = json.load(open(p))
        for r in d.get("models", d.get("records", [])) if isinstance(d, dict) else d:
            r = dict(r)
            r.setdefault("model", d.get("model") if isinstance(d, dict) else None)
            P, G, T = num(r.get("power_kW")), num(r.get("exh_flow_kg_s")), num(r.get("exh_temp_C"))
            if not (P and G and T):
                continue
            var = str(r.get("variant", "")).upper()
            gen = var.startswith(("GS", "PG"))
            Psh = P / ETA_GEN if gen else P
            out.append(dict(model=r["model"], variant="генератор" if gen else "привод", P_kW=P, P_shaft_kW=round(Psh),
                            G_air=round(G / (1 + FAR), 2), EGT0=T, eta0=num(r.get("eff_pct")), rpm0=rpm_of(r),
                            shafts0=r.get("shafts"), PR0=num(r.get("PR")), comp0=r.get("comp_stages"), turb0=r.get("turb_stages"),
                            pkg=(num(r.get("pkg_L_m")), num(r.get("pkg_W_m")), num(r.get("pkg_H_m"))),
                            w0=round(Psh / (G / (1 + FAR)), 1), src=os.path.basename(p),
                            wheels0=stages_sum(r.get("comp_stages")) + stages_sum(r.get("turb_stages")),
                            recup=any(k in str(r["model"]) for k in RECUP)))
    return out


def best_for(slot, pts, ak, ck, alpha, dT=None):
    Gcap = slot["G_air"] * (1 + alpha)
    w_req = slot["P_shaft_kW"] / Gcap
    cand = [p for p in pts if p["arch"] == ak and p["cool"] == ck and p["w"] >= w_req
            and (dT is None or p["EGT"] <= slot["EGT0"] + dT)]
    if not cand:
        return None
    Tmin = min(p["T41"] for p in cand)
    p = max((p for p in cand if p["T41"] == Tmin), key=lambda p: p["eta"])
    G = slot["P_shaft_kW"] / p["w"]
    res = dict(arch=ak, cool=ck, PR=p["PR"], T41=p["T41"], w=round(p["w"], 1), eta=round(p["eta"], 1), G=round(G, 2),
               dG=round((G / slot["G_air"] - 1) * 100, 1), EGT=round(p["EGT"]), dEGT=round(p["EGT"] - slot["EGT0"]),
               fuel_pct=round((slot["eta0"] / p["eta"] - 1) * 100) if slot["eta0"] else None)
    # колесо компрессора (первое): диаметр и собственная частота
    PRs1 = p["PR"] if ARCH[ak]["nc"] == 1 else p["PR"] ** 0.5
    st = E.comp_stage(E.T0, E.P0, PRs1, G, "ti")
    res.update(D2=round(st["D2"], 3), n_gg=round(st["rpm"]), U2=round(st["U2"]))
    if st["D2"] <= D2_MAX:
        res["wheel"] = "односторонний вход"
    elif st["D2"] / math.sqrt(2) <= D2_MAX:
        res.update(wheel="двусторонний вход", D2=round(st["D2"] / math.sqrt(2), 3), n_gg=round(st["rpm"] * math.sqrt(2)))
    else:
        ngg = math.ceil((st["D2"] / math.sqrt(2) / D2_MAX) ** 2)
        res.update(wheel=f"{ngg} газогенератора параллельно, двусторонний вход", D2=round(st["D2"] / math.sqrt(2 * ngg), 3),
                   n_gg=round(st["rpm"] * math.sqrt(2 * ngg)), n_gg_units=ngg)
        if ak in ("1В-1Ц", "1В-2Ц"):
            return None
    # габарит: корпус центробежного компрессора ~1,9 D2 плюс 0,5 м на обвязку против меньшего из ширины и высоты пакета
    W, H = slot["pkg"][1], slot["pkg"][2]
    env = 1.9 * res["D2"] + 0.5
    res["envelope_m"] = round(env, 2)
    res["envelope_ok"] = None if not (W and H) else bool(env <= min(W, H))
    if ak in ("1В-1Ц", "1В-2Ц"):
        # частота вала одновального двигателя должна совпасть с частотой входа редуктора пакета
        need = slot["rpm0"]
        if not need:
            return None
        k = need / st["rpm"]
        # коэффициент расхода 0,06...0,11 вместо 0,085 -> частота x0,84...1,14; двусторонний вход -> x1,41
        if 0.84 <= k <= 1.14:
            res.update(rpm_fit="прямо", rpm_k=round(k, 2), eta=round(p["eta"] - 1.5 * abs(math.log(k)) / math.log(1.14), 1))
        elif 0.84 * 1.414 <= k <= 1.14 * 1.414:
            if res["wheel"] == "односторонний вход":
                res.update(wheel="двусторонний вход", D2=round(st["D2"] / math.sqrt(2), 3), n_gg=round(st["rpm"] * math.sqrt(2)))
                res["envelope_m"] = round(1.9 * res["D2"] + 0.5, 2)
                res["envelope_ok"] = None if not (W and H) else bool(res["envelope_m"] <= min(W, H))
            res.update(rpm_fit="двусторонний вход колеса", rpm_k=round(k, 2),
                       eta=round(p["eta"] - 0.5 - 1.5 * abs(math.log(k / 1.414)) / math.log(1.14), 1))
        else:
            return None
        res["n_out"] = round(need)
    else:
        # силовая турбина под исходную частоту: окружная скорость <= 380 м/с, нагрузка ступени <= 1,6
        need = slot["rpm0"] or 9500
        dh = p["w"] * 1e3 / 0.99 * (1.0)          # работа СТ на кг воздуха ~ удельная работа цикла
        U = None
        for z in (1, 2, 3):
            Uz = math.sqrt(dh / z / 1.6)
            if Uz <= 380:
                U, zpt = Uz, z
                break
        if U is None:
            U, zpt = 380.0, 3
        res.update(pt_stages=zpt, pt_D=round(60 * U / (math.pi * need), 2), n_out=round(need))
    res["parts"] = dict(rotors=ARCH[ak]["rotors"], wheels=ARCH[ak]["nc"] + (1 if ARCH[ak]["rotors"] == 1 else 1 + res.get("pt_stages", 0)),
                        cooled=COOL[ck]["rank"])
    return res


def choose(slot, pts, alpha, dT=None):
    order = []
    for ak in (["1В-1Ц", "1В-2Ц"] if slot["variant"] == "генератор" else []) + ["2В-1Ц", "2В-2Ц"]:
        for ck in ("U", "N", "K"):
            order.append((ak, ck))
    # сортировка по сложности: охлаждение важнее числа роторов для надёжности; затем роторы; затем колёса
    order.sort(key=lambda t: (COOL[t[1]]["rank"], ARCH[t[0]]["rotors"], ARCH[t[0]]["nc"]))
    for ak, ck in order:
        r = best_for(slot, pts, ak, ck, alpha, dT)
        if r:
            return r
    return None


def frames(chosen, span=0.75):
    """минимальное число геометрических типоразмеров газогенератора: типоразмер с расходом Gf обслуживает
    установки с расходом от span*Gf до Gf (дросселирование температурой газа)"""
    Gs = sorted(c["G"] for c in chosen if c)
    out, i = [], len(Gs) - 1
    while i >= 0:
        Gf = Gs[i]
        members = [g for g in Gs if span * Gf <= g <= Gf]
        out.append(dict(G_frame=round(Gf, 1), n=len(members), G_min=round(min(members), 1)))
        i -= len(members)
        while i >= 0 and Gs[i] >= span * Gf:
            i -= 1
    return out[::-1]


if __name__ == "__main__":
    paths = sys.argv[1:] or sorted(glob.glob("fleet/*.json"))
    S = slots(paths)
    pts = grid()
    res = []
    for s in S:
        row = dict(s)
        for k, sc in SCEN.items():
            row[k] = None if s["recup"] else choose(s, pts, sc["alpha"], sc["dT"])
        row["без_ограничения_выхлопа"] = None if s["recup"] else choose(s, pts, 0.10, None)
        res.append(row)
    fam = {k: frames([r[k] for r in res]) for k in list(SCEN) + ["без_ограничения_выхлопа"]}
    fam_deep = {k: frames([r[k] for r in res], span=0.6) for k in SCEN}
    meta = {"method": "0D-цикл (gtu_cycle.cycle, Cantera GRI-Mech 3.0, метан, ISO) + одномерный размер колеса (gtu_explore.comp_stage); "
                      "отбор схем по минимуму движущихся частей и охлаждаемых деталей",
            "software": "cantera " + __import__("cantera").__version__ + ", python3",
            "date": __import__("datetime").date.today().isoformat(), "inputs": [os.path.relpath(p) for p in paths]}
    json.dump(dict(**meta, arch=ARCH, cool=COOL, scenarios=SCEN, d2_max=D2_MAX, slots=res, frames=fam, frames_span06=fam_deep), open("gtu_family_results.json", "w"),
              ensure_ascii=False, indent=1, default=float)
    for r in res:
        line = f"{r['model'][:22]:22s} {r['variant'][:3]} P={r['P_kW']:>6.0f} G={r['G_air']:>5.1f} w0={r['w0']:>5.0f} n={r['rpm0']} |"
        for k in list(SCEN) + ["без_ограничения_выхлопа"]:
            c = r[k]
            line += (f" {k[:4]}: {c['arch']}/{c['cool']} T{c['T41']} π{c['PR']} η{c['eta']} G{c['dG']:+.0f}% ΔT{c['dEGT']:+d} D{c['D2']} |" if c else f" {k[:4]}: —  |")
        print(line)
    print(json.dumps(fam, ensure_ascii=False))
