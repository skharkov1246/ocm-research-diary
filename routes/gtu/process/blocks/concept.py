"""Описание концепции «от топлива до вала»: шесть блоков + сборка, значения по умолчанию, перевод в концепцию воронки.

Концепция - словарь по блокам (пропущенное берётся из DEFAULT):
  fuel        - подготовка топлива: prep = нет | подогрев | паровой риформинг | сухой риформинг;
  compressor  - arch 1Ц | 2Ц | Д | ОЦ (+ n_ax), inlet одно | двух, mat Ti64 | Al | Ni, make поковка | литьё | печать, pi, split,
                mu, промежуточное охлаждение intercool нет | воздух (после ступени ic_after);
  combustor   - scheme LPP | FLOX | TVC | CAT | RDC | ВР, chamber кольцевая | выносная, dp_source comb_l0 | track | comb_l0_min;
  turbine     - type О | Р, nt (auto), cool U | N | K | P | C, mat, make, hot_make, hot_mat, t41, npt (auto) силовой турбины;
  hx          - recup (0 или ε 0,70-0,92), recup_make фольга | печать, recup_integral, reheat, x_rh, t_rh, ibc;
  shaft       - drive прямой | свободная СТ | быстроходный (через существующий редуктор пакета, только генераторный);
  build       - rotor_join, bear масло | магнит, case картридж | оболочка | разъём, conn, cores, comb_make, case_make, kD_case;
  eta         - source воронка | задано (comp [...], gg, pt) - КПД узлов: из моделей воронки или заданные (калибровка, сверка);
  options     - rec_model enthalpy | temperature, T_rec_max (газ на входе в рекуператор), dp_inlet.
"""
import copy
import math

from . import _track as T

DEFAULT = dict(
    kv="KV-3", name=None,
    fuel=dict(prep="нет", T_max=None, sc=2.0, R_egr=1.0),
    compressor=dict(arch="2Ц", n_ax=0, inlet="одно", mat="Ti64", make="поковка", pi=None, split=0.5, mu=None, kcm=0.28, zc=11,
                    intercool="нет", ic_after=None),
    combustor=dict(scheme="LPP", chamber="кольцевая", dp_source="comb_l0", pgc_real=None, dp=None),
    turbine=dict(type="О", nt="auto", cool="U", mat="SC", make="литьё", hot_make="литьё", hot_mat="K438", t41=1000.0, lt=None,
                 phit=0.5, npt="auto", kUpt=1.0, phipt=0.5, zt=None, zn=None),
    hx=dict(recup=0.0, recup_make="фольга", recup_integral=0, reheat=False, x_rh=0.5, t_rh=None, ibc=False),
    shaft=dict(drive="прямой", phi1=None, n_out=None),
    build=dict(rotor_join="сборный", bear="масло", case="картридж", conn="болты", cores=1, comb_make="обработка", case_make="литьё",
               split_bl=0, kD_case=1.9),
    eta=dict(source="воронка"),
    options=dict(rec_model="enthalpy", T_rec_max=None),
)

PI_DEFAULT = {"1Ц": 6.0, "Д": 4.5, "2Ц": 13.5, "ОЦ": 12.0}          # середины диапазонов воронки (funnel_space.PI_RANGE) - оценка
MU_DEFAULT = {"Ц": 0.70, "Д": 0.60, "О": 0.28}                      # comp_stage (gtu_explore) - 0,70; середины MU_RANGE воронки
LT_DEFAULT = {"О": 1.2, "Р": 0.70}                                  # ψ при наибольшем КПД eta_axial; U/C0 оптимума turb_radial
PSI_AUTO = 1.6           # число ступеней осевой турбины «auto»: наименьшее с ψ <= 1,6 при пределе окружной скорости (оценка по
#                          frame3: при ψ 1,79 число Маха за сопловым аппаратом 1,24-1,44 > 1,2 - frame3_point_results.json turbine)
NU_AUTO = 0.60           # радиальная: наименьшее число ступеней с U/C0 >= 0,60 при пределе окружной скорости (оценка)
ZT = {"О": 60, "Р": 13}  # лопаток колеса турбины (воронка: осевая 30-80, радиальная 9-17) - только для счёта деталей, оценка
ZN = {"О": 40, "Р": 18}  # лопаток соплового аппарата (воронка: осевая 20-60, радиальная 12-24) - оценка
SHAFT_DRIVES = ("прямой", "свободная СТ", "быстроходный")
PREP = ("нет", "подогрев", "паровой риформинг", "сухой риформинг")
SCHEMES = ("LPP", "FLOX", "TVC", "CAT", "RDC", "ВР")


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k].update(v)
        else:
            out[k] = v
    return out


def stage_types(C):
    c = C["compressor"]
    return {"1Ц": ["Ц"], "Д": ["Д"], "2Ц": ["Ц", "Ц"], "ОЦ": ["О"] * max(c["n_ax"], 1) + ["Ц"]}[c["arch"]]


def normalize(concept):
    """концепция с заполненными значениями по умолчанию; ошибки схемы - ValueError"""
    C = _merge(DEFAULT, concept)
    c, cb, t, hx, sh, b = C["compressor"], C["combustor"], C["turbine"], C["hx"], C["shaft"], C["build"]
    if c["arch"] not in PI_DEFAULT:
        raise ValueError(f"компрессор: архитектура {c['arch']}")
    if c["arch"] == "ОЦ" and c["n_ax"] < 1:
        c["n_ax"] = 1
    if c["arch"] != "ОЦ":
        c["n_ax"] = 0
    if c["pi"] is None:
        c["pi"] = PI_DEFAULT[c["arch"]]
    typ = stage_types(C)
    if c["mu"] is None:
        c["mu"] = [MU_DEFAULT[x] for x in typ]
    elif isinstance(c["mu"], (int, float)):
        c["mu"] = [float(c["mu"])] * len(typ)
    if len(c["mu"]) != len(typ):
        raise ValueError(f"компрессор: mu на {len(typ)} ступени")
    if c["intercool"] not in ("нет", "воздух"):
        raise ValueError(f"промохлаждение: {c['intercool']}")
    if c["ic_after"] is None:
        c["ic_after"] = len(typ) - 2 if len(typ) >= 2 else -1
    if cb["scheme"] not in SCHEMES:
        raise ValueError(f"камера: схема {cb['scheme']}")
    if C["fuel"]["prep"] not in PREP:
        raise ValueError(f"топливо: {C['fuel']['prep']}")
    if sh["drive"] not in SHAFT_DRIVES:
        raise ValueError(f"вал: {sh['drive']}")
    if t["lt"] is None:
        t["lt"] = LT_DEFAULT[t["type"]]
    if t["zt"] is None:
        t["zt"] = ZT[t["type"]]
    if t["zn"] is None:
        t["zn"] = ZN[t["type"]]
    if hx["t_rh"] is None:
        hx["t_rh"] = t["t41"]
    if hx["recup"] and not (0.5 <= hx["recup"] < 1.0):
        raise ValueError(f"рекуператор: ε {hx['recup']}")
    if b["case"] not in ("картридж", "оболочка", "разъём"):
        raise ValueError(f"корпус: {b['case']}")
    C["rotors"] = 2 if sh["drive"] == "свободная СТ" else 1
    return C


def label(C):
    """короткое имя концепции из блоков"""
    if C.get("name"):
        return C["name"]
    c, t, hx, cb = C["compressor"], C["turbine"], C["hx"], C["combustor"]
    s = [f"{C['shaft']['drive']}", f"{c['arch']}{c['n_ax'] if c['arch'] == 'ОЦ' else ''} π{c['pi']:.3g}",
         f"{'осевая' if t['type'] == 'О' else 'радиальная'} ×{t.get('nt_used') or t['nt']} {t['cool']}", f"T41 {t['t41']:.0f}", cb["scheme"]]
    if c["intercool"] != "нет":
        s.append("промохлаждение")
    if hx["recup"]:
        s.append(f"рекуператор {hx['recup']:.2f} ({hx['recup_make']})")
    if hx["reheat"]:
        s.append("подогрев")
    if hx["ibc"]:
        s.append("IBC")
    if C["fuel"]["prep"] != "нет":
        s.append(C["fuel"]["prep"])
    return ", ".join(s)


def funnel_arch(C):
    """архетип воронки - только метка (ниша); на модели деталей, обслуживания, Plug&Play и изготовимости не влияет"""
    t, b, cb, hx = C["turbine"], C["build"], C["combustor"], C["hx"]
    if cb["scheme"] == "RDC":
        return "RDC"
    if cb["scheme"] == "ВР":
        return "волновой ротор"
    if b["cores"] > 1:
        return "кластер"
    if t["cool"] == "C":
        return "керамика"
    if hx["recup"]:
        return "печатный рекуператор"
    if b["bear"] == "магнит":
        return "безмасляный"
    if b["rotor_join"] == "заодно" and t["make"] == "печать":
        return "печатный монолит"
    if b["rotor_join"] in ("заодно", "сварка") and t["make"] == "литьё":
        return "литой монолит"
    if t["hot_make"] == "печать":
        return "гибрид"
    if b["case"] == "разъём":
        return "традиционная"
    return "гибрид"


def to_funnel(C, tgt, nt=None, npt=None, pgain=0.0):
    """концепция воронки (все оси funnel_space.CAT и CONT) из концепции процесса. Ось «корпус» - funnel_space.apply_case."""
    FS = T.FS()
    c, cb, t, hx, sh, b = C["compressor"], C["combustor"], C["turbine"], C["hx"], C["shaft"], C["build"]
    typ = stage_types(C)
    mus = c["mu"]
    if c["arch"] == "ОЦ":
        mu1, mu2 = mus[0], mus[-1]
    else:
        mu1, mu2 = mus[0], (mus[1] if len(mus) > 1 else mus[0])
    k0 = next(k for k, x in enumerate(typ) if x != "О")
    comb = {"RDC": "RDC", "ВР": "ВР"}.get(cb["scheme"], cb["chamber"])
    fc = dict(idx=-1, n=0, seed=FS.SEED, kv=tgt["kv"], arch=funnel_arch(C), rotors=C["rotors"], comp=c["arch"], n_ax=c["n_ax"],
              inlet=c["inlet"], turb=t["type"], nt=int(nt or (t["nt"] if t["nt"] != "auto" else 1)),
              npt=int(npt or (t["npt"] if t["npt"] != "auto" else 1)) if C["rotors"] == 2 else 0,
              cool=t["cool"], recup="печатный" if hx["recup"] else "нет", comb=comb, bear=b["bear"], cores=b["cores"],
              rotor_join=b["rotor_join"], m_cw=c["make"], m_tw=t["make"], m_hot=t["hot_make"], m_comb=b["comb_make"],
              m_case=b["case_make"], mat_cw=c["mat"], mat_tw=t["mat"], mat_hot=t["hot_mat"], hot_integral=0,
              recup_integral=int(hx["recup_integral"]), conn=b["conn"], cartridge=1 if b["case"] == "картридж" else 0,
              split_bl=b["split_bl"], mu1=mu1, mu2=mu2, lt=t["lt"], eps_rec=float(hx["recup"] or 0.0), pgain=float(pgain),
              kcm=c["kcm"], zc=c["zc"], zt=t["zt"], zn=t["zn"], split=c["split"], t41=float(t["t41"]), pi=float(c["pi"]),
              phi1=float(sh["phi1"] or FS.PHI_OPT["Д" if typ[k0] == "Д" else "Ц"]) if C["rotors"] == 2 else 0.0,
              phit=t["phit"], kUpt=t["kUpt"] if C["rotors"] == 2 else 0.0, phipt=t["phipt"] if C["rotors"] == 2 else 0.0,
              kD_case=b["kD_case"])
    case = b["case"]
    return FS.apply_case(fc, case)


def stages_pi(C):
    """степени сжатия ступеней (как funnel_space.stages_pi: деление ln π долей split)"""
    c = C["compressor"]
    pi, s = c["pi"], c["split"]
    if c["arch"] in ("1Ц", "Д"):
        return [pi]
    if c["arch"] == "2Ц":
        return [pi ** s, pi ** (1 - s)]
    na = c["n_ax"]
    return [pi ** (s / na)] * na + [pi ** (1 - s)]


def split_for_pi1(pi, pi1):
    """доля ln π первой ступени по её степени сжатия"""
    return math.log(pi1) / math.log(pi)
