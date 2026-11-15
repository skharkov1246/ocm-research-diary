"""Извлечение обобщённых характеристик компрессора из пакета pyCycle (NASA Glenn, Apache-2.0).

Зачем. Для согласования узлов ГТД вне ISO нужна характеристика компрессора двигателя, которой
у нас нет (закрытые данные изготовителя). Берём открытую обобщённую характеристику: карты NPSS,
переведённые в pyCycle 4.4.0 (`pycycle/maps/HPC_map.py` — «HBTF HPC map from NPSS»,
`pycycle/maps/axi5.py` — «axi-5.map Compressor map from NPSS»), и масштабируем на расчётную точку
Taurus 60 (gtu_shaft.CompressorMap). Это ДОПУЩЕНИЕ: форма линий частоты чужого компрессора.

Как повторить (из корня репозитория):
    d=$(mktemp -d) && pip download --no-deps om-pycycle==4.4.0 -d "$d"
    python3 -I routes/gtu/shaft/extract_pycycle_maps.py "$d"/om_pycycle-4.4.0-py3-none-any.whl
Скрипт не импортирует и не выполняет код пакета: файл карты читается из архива как текст и
разбирается `ast.literal_eval` (только числовые литералы). Результат — maps/<имя>.json рядом.
"""
import ast, hashlib, json, os, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
MAPS = {   # имя в репозитории: (путь в архиве, переменная карты)
    "npss_hbtf_hpc": ("pycycle/maps/HPC_map.py", "HPCMap"),
    "npss_axi5": ("pycycle/maps/axi5.py", "AXI5"),
}
ARRAYS = ("NcMap", "RlineMap", "alphaMap", "WcMap", "effMap", "PRmap")


def _literal(node):
    """np.array([...]) -> список; иначе литерал"""
    if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "array" and node.args:
        return ast.literal_eval(node.args[0])
    return ast.literal_eval(node)


def parse(src, var):
    tree = ast.parse(src)
    out, defaults, doc = {}, {}, None
    for st in tree.body:
        if isinstance(st, ast.Expr) and isinstance(getattr(st, "value", None), ast.Constant) \
                and isinstance(st.value.value, str) and doc is None:
            doc = st.value.value.strip()
        if not isinstance(st, ast.Assign) or len(st.targets) != 1:
            continue
        t = st.targets[0]
        if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name) and t.value.id == var \
                and t.attr in ARRAYS + ("RlineStall",):
            out[t.attr] = _literal(st.value)
        if isinstance(t, ast.Subscript) and isinstance(t.value, ast.Attribute) and t.value.attr == "defaults":
            defaults[ast.literal_eval(t.slice)] = ast.literal_eval(st.value)
    return out, defaults, doc


def main(wheel):
    wsha = hashlib.sha256(open(wheel, "rb").read()).hexdigest()
    with zipfile.ZipFile(wheel) as z:
        lic = z.read("om_pycycle-4.4.0.dist-info/licenses/LICENSE.txt").decode("utf-8", "replace")
        for name, (path, var) in MAPS.items():
            raw = z.read(path)
            arr, defaults, doc = parse(raw.decode("utf-8"), var)
            missing = [k for k in ARRAYS if k not in arr]
            if missing:
                raise SystemExit(f"{path}: нет массивов {missing}")
            rec = dict(
                name=name, description=doc, map_variable=var,
                source=dict(package="om-pycycle 4.4.0 (PyPI), NASA Glenn Research Center, pyCycle",
                            wheel=os.path.basename(wheel), wheel_sha256=wsha, file=path,
                            file_sha256=hashlib.sha256(raw).hexdigest(),
                            project="https://github.com/OpenMDAO/pyCycle",
                            paper="Hendricks E.S., Gray J.S. pyCycle: A Tool for Efficient Optimization of Gas "
                                  "Turbine Engine Cycles. Aerospace 2019, 6(8), 87 (ссылка из открытых сведений о "
                                  "пакете; статья из контейнера не открывалась)"),
                license="Apache-2.0 (текст лицензии пакета: '" + " ".join(lic.split())[:160] + "...')",
                axes=dict(alphaMap="угол поворотных лопаток карты; используется срез alpha=0",
                          NcMap="приведённая частота, отн. (в pyCycle подписано 'rpm', значения относительные)",
                          RlineMap="линия R: 1 — граница помпажа (RlineStall), 3 — сторона запирания"),
                units=dict(WcMap="lbm/s (после масштабирования не важно)", effMap="изоэнтропический КПД",
                           PRmap="степень повышения давления"),
                defaults=defaults, RlineStall=arr.get("RlineStall"),
                **{k: arr[k] for k in ARRAYS})
            dst = os.path.join(HERE, "maps", name + ".json")
            json.dump(rec, open(dst, "w"), ensure_ascii=False, indent=1)
            print(dst, "Nc", len(arr["NcMap"]), "R", len(arr["RlineMap"]))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
