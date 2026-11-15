"""Доступ к модулям трека «Простая литая ГТУ-заменитель» без их правки.

Модули импортируются как есть: gtu_cycle и gtu_explore при импорте делают chdir в папку трека, поэтому рабочая папка
возвращается обратно. Ничего в чужих файлах не меняется; единственное вмешательство — регистрация расчётной точки
в кэше funnel_space._TGT (словарь точек модуля воронки) для типоразмеров и частот вала, которых нет в funnel_space.KV
(функции воронки берут точку по ключу concept['kv']). Ключи регистрируются с префиксом PB: и не совпадают с KV-1/3/6.
"""
import importlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROCESS = os.path.dirname(HERE)
GTU = os.path.dirname(PROCESS)
REPO = os.path.abspath(os.path.join(GTU, "..", ".."))
FUNNEL = os.path.join(GTU, "funnel")
CYCLES = os.path.join(GTU, "organic", "cycles")
COMB = os.path.join(GTU, "combustor")
FLEET = os.path.join(GTU, "fleet")

F_FAMILY = os.path.join(GTU, "gtu_family_results.json")
F_FRAME3 = os.path.join(GTU, "frame3", "frame3_point_results.json")
F_CYCLES = os.path.join(CYCLES, "cycles_results.json")
F_COMB_TARGETS = os.path.join(COMB, "targets.json")
F_COMB_L0 = os.path.join(COMB, "l0_aws", "comb_l0_KV-3_n262144_results.json")

_MOD = {}


def rel(p):
    return os.path.relpath(p, REPO)


def _imp(name, path):
    if name in _MOD:
        return _MOD[name]
    cwd = os.getcwd()
    if path not in sys.path:
        sys.path.insert(0, path)
    try:
        m = importlib.import_module(name)
    finally:
        os.chdir(cwd)
    _MOD[name] = m
    return m


def FS():
    """воронка, уровень 0: пространство, быстрый цикл, геометрия, правила (funnel/funnel_space.py)"""
    return _imp("funnel_space", FUNNEL)


def L0():
    """воронка, уровень 0: решение с поправками КПД, запись, оценка (funnel/funnel_l0.py)"""
    return _imp("funnel_l0", FUNNEL)


def MK():
    return _imp("make_model", FUNNEL)


def PM():
    return _imp("parts_model", FUNNEL)


def SM():
    return _imp("service_model", FUNNEL)


def PP():
    return _imp("plugplay_model", FUNNEL)


def CY():
    """перебор циклов (organic/cycles/cycles.py): Cantera, 576 топологий, теплообменники, TRL, детали"""
    return _imp("cycles", CYCLES)


def CS():
    """воронка камеры сгорания (combustor/comb_space.py): схемы, TRL схем"""
    return _imp("comb_space", COMB)


def GC():
    """модель цикла трека (gtu_cycle.py: cycle, cycle_r, recup)"""
    return _imp("gtu_cycle", GTU)


def E():
    """свойства газа: таблицы gtu_explore (Cantera GRI-Mech 3.0) через funnel_space.FastGas"""
    return FS()._E()


_JSON = {}


def load_json(path):
    if path not in _JSON:
        with open(path, encoding="utf-8") as f:
            _JSON[path] = json.load(f)
    return _JSON[path]


# =========================================================================== точка пакета (стык)
FLEET_FILE = {"Taurus 60": "gtu_fleet_taurus60.json", "Saturn 20": "gtu_fleet_saturn20.json", "Titan 130": "gtu_fleet_titan130.json",
              "Centaur 40": "gtu_fleet_centaur40.json", "Centaur 50": "gtu_fleet_centaur50.json", "Taurus 65": "gtu_fleet_taurus65.json",
              "Taurus 70": "gtu_fleet_taurus70.json", "Mars 90": "gtu_fleet_mars90.json", "Mars 100": "gtu_fleet_mars100.json",
              "Mercury 50": "gtu_fleet_mercury50.json", "SGT-100": "gtu_fleet_sgt100.json", "SGT-300": "gtu_fleet_sgt300.json",
              "SGT-400": "gtu_fleet_sgt400.json", "SGT-500": "gtu_fleet_sgt500.json", "SGT-A05": "gtu_fleet_sgta05.json"}


def package_target(kv="KV-3", n_out=None):
    """точка пакета: для KV-1/3/6 - funnel_space.target (gtu_family_results.json + паспорт fleet) без изменений;
    n_out - частота вала двигателя, отличная от частоты пакета (быстроходный вариант через редуктор пакета): копия точки
    регистрируется в кэше воронки под своим ключом, n_pkg_true хранит частоту пакета. Плюс расход выхлопа пакета
    (паспорт fleet, исполнение генераторное) - для допуска +20 % по выхлопу."""
    fs = FS()
    t0 = fs.target(kv)
    t = dict(t0)
    t["n_pkg_true"] = t0["n_pkg"]
    t["W_exh0"] = _exh_flow(t0["model"])
    if n_out is not None and abs(n_out - t0["n_pkg"]) > 0.5:
        key = f"PB:{kv}:n{int(round(n_out))}"
        t.update(kv=key, n_pkg=float(n_out))
        fs._TGT[key] = t
        return t
    t["kv"] = kv
    return t


def _exh_flow(model):
    fn = FLEET_FILE.get(model)
    if not fn:
        return None
    d = load_json(os.path.join(FLEET, fn))
    gen = [r for r in d["records"] if r.get("variant") in ("GS", "PG") or str(r.get("variant", "")).startswith("PG")]
    rec = gen[0] if gen else d["records"][0]
    v = rec.get("exh_flow_kg_s")
    return float(v) if v else None
