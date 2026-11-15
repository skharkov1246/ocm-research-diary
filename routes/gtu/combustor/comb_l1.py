"""Воронка камеры сгорания KV, уровень 1: одна часть - осесимметричный реагирующий RANS (OpenFOAM) по вариантам.

Вход - l1_inputs_<kv>.json (prepare_l1.py). Часть shard из shards получает варианты с номером i (позиция во входном
файле, i < n) при i % shards == shard - разбиение одно, как у уровня 0. На вариант четыре кейса (comb_cfd.py):
ячейка на 100, 50, 30 % и труба на 100 % (вход трубы - выход ячейки 100 %); у каждого k-го варианта ещё ячейка 100 %
на сетке medium (проверка сетки). Кейс - одно ядро; кейсы части идут пулом на --workers процессов (на c7a.16xlarge -
$NPROC = 64). Итоги кейса - comb_post.analyse; в папку части - JSON части и сжатые итоги кейсов (логи, case.json,
последняя запись основных полей), без полного набора полей.

Запуск:
  python comb_l1.py --inputs l1_inputs_KV-3.json --n 500 --shard 3 --shards 25 --out "$GTU_OUT" --workers "$NPROC"
  python comb_l1.py --inputs /tmp/l1.json --n 2 --shard 0 --shards 1 --out /tmp/l1 --of 1912 --iters 200,100,100
  python comb_l1.py ... --generate-only      # только кейсы и сетки (blockMesh, checkMesh), без решателя
AWS не вызывается; нужен OpenFOAM в окружении (reactingFoam, blockMesh, setFields, foamDictionary).
"""
import argparse, hashlib, json, os, platform, shutil, subprocess, sys, tarfile, time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_cfd as CF   # noqa: E402
import comb_post as CP  # noqa: E402

SCHEMA = "kv.combustor.l1.shard/1"
KEEP_FIELDS = ("T", "U", "p", "CH4", "O2", "CO", "CO2", "H2O", "NO", "NO2", "OH", "k", "epsilon")


def shard_indices(n, shard, shards):
    if not (shards >= 1 and 0 <= shard < shards and n >= 1):
        raise ValueError(f"неверные n={n}, shard={shard}, shards={shards}")
    return list(range(shard, n, shards))


def code_hash():
    h = hashlib.sha256()
    for f in ("comb_cfd.py", "comb_post.py", "comb_l1.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    for f in ("reactionsGRI", "thermo.compressibleGasGRI", "radiationProperties"):
        h.update(open(os.path.join(CF.FOAM_RES, f), "rb").read())
    return h.hexdigest()[:16]


def shard_name(kv, n, shard, shards):
    return f"comb_l1_{kv}_n{n}_s{shard:03d}of{shards:03d}.json"


def _run(cmd, cwd, log, timeout_s):
    t0 = time.time()
    with open(os.path.join(cwd, log), "w") as f:
        try:
            r = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, timeout=timeout_s)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = 124
    return rc, time.time() - t0


def run_case(case, hot_init, iters, timeout_s, generate_only=False):
    """тот же порядок, что Allrun кейса: сетка -> без химии -> горячее поле -> химия; коды и время шагов"""
    steps = []

    def step(cmd, log, frac=1.0):
        left = max(60.0, timeout_s - sum(s["s"] for s in steps))
        rc, dt = _run(cmd, case, log, left * frac)
        steps.append(dict(cmd=" ".join(cmd), rc=rc, s=round(dt, 2)))
        return rc
    if step(["blockMesh"], "log.blockMesh") != 0:
        return steps, "сетка не построена"
    step(["checkMesh"], "log.checkMesh")
    mesh_ok = "Mesh OK." in open(os.path.join(case, "log.checkMesh"), errors="replace").read()
    if not mesh_ok:
        return steps, "checkMesh: сетка с ошибками"
    if generate_only:
        return steps, "только сетка"
    fd = lambda d, e, v: step(["foamDictionary", d, "-entry", e, "-set", str(v)], "log.foamDictionary")
    fd("constant/chemistryProperties", "chemistry", "off")
    fd("system/controlDict", "endTime", iters["cold"])
    if step(["reactingFoam"], "log.reactingFoam.cold", 0.4) != 0:
        return steps, "решатель упал (без химии)"
    if hot_init:
        step(["setFields"], "log.setFields")
    fd("constant/chemistryProperties", "chemistry", "on")
    fd("system/controlDict", "endTime", iters["cold"] + iters["hot"])
    rc = step(["reactingFoam"], "log.reactingFoam.hot")
    if rc == 124:
        return steps, "не хватило времени (сторож кейса)"
    if rc != 0:
        return steps, "решатель упал (с химией)"
    return steps, "ok"


def _pack(case, dest):
    """сжатые итоги кейса: логи, case.json, result.json, система, основная часть последней записи"""
    ts = CP.latest_times(case, 1)
    with tarfile.open(dest, "w:gz") as tz:
        for f in os.listdir(case):
            p = os.path.join(case, f)
            if f.startswith("log.") or f in ("case.json", "result.json", "Allrun"):
                tz.add(p, arcname=f)
        tz.add(os.path.join(case, "system"), arcname="system")
        pm = os.path.join(case, "constant", "polyMesh")
        if os.path.isdir(pm):
            tz.add(pm, arcname="constant/polyMesh")
        if ts:
            for fld in KEEP_FIELDS:
                p = os.path.join(case, ts[0], fld)
                if os.path.exists(p):
                    tz.add(p, arcname=f"{ts[0]}/{fld}")


def task(args):
    """один кейс: (вариант, вид, нагрузка, сетка) -> итог"""
    e, kind, li, level, work, of, iters, timeout_s, gen_only, cell_out = args
    t0 = time.time()
    name = f"{kind}_{int(round(100 * e['loads'][li]['load']))}_{level}"
    case = os.path.join(work, f"v{e['_pos']:04d}_i{e['idx']}", name)
    res = dict(name=name, kind=kind, load=e["loads"][li]["load"], level=level)
    try:
        if kind == "cell":
            desc = CF.build_cell_case(e, li, case, level, of, iters)
            hot = True
        else:
            desc = CF.build_liner_case(e, cell_out, case, level, of, iters)
            hot = False
        res["n_cells"] = desc["n_cells"]
        steps, st = run_case(case, hot, iters, timeout_s, gen_only)
        res.update(steps=steps, run_status=st)
        if st == "ok":
            a = CP.analyse(case, dict(T3=e["loads"][li]["T3_K"], lhv=e["lhv"]))
            json.dump(a, open(os.path.join(case, "result.json"), "w"), ensure_ascii=False, default=float)
            res["result"] = a
        elif st == "только сетка":
            res["result"] = dict(n_cells=desc["n_cells"], patches=desc["patches"])
    except Exception as ex:          # кейс не теряется: статус «ошибка» с причиной
        res.update(run_status="ошибка", error=f"{type(ex).__name__}: {str(ex)[:300]}")
    res["wall_s"] = round(time.time() - t0, 1)
    return e["_pos"], res, case


def cell_outlet_state(cell_res, e):
    """осреднённое состояние выхода ячейки 100 % для входа трубы (на одну трубу)"""
    a = cell_res["result"]
    ms = a["outlet_mixed_state"]
    c, L = e["concept"], e["loads"][0]
    n_el_can = c["n_el"] / e["n_can"]
    m_el = a["mass_out"] * (3.141592653589793 / CF.TAN_A)            # клин -> полный элемент
    _ = c, L
    return dict(T=ms["T"], Y=ms["Y"], m_dot_can=m_el * n_el_can)


def run(inputs, n, shard, shards, out, workers, level="coarse", of="2512", iters=None, grid_every=10,
        timeout_min=180, generate_only=False, keep_work=False):
    t0 = time.time()
    data = json.load(open(inputs, encoding="utf-8"))
    ents = data["entries"]
    if n > len(ents):
        raise SystemExit(f"во входном файле {len(ents)} вариантов, запрошено {n}")
    it = dict(CF.ITERS, **(iters or {}))
    idx = shard_indices(n, shard, shards)
    work = os.path.join(out, f"work_s{shard:03d}")
    os.makedirs(work, exist_ok=True)
    for i in idx:
        ents[i]["_pos"] = i
    tasks = []
    for i in idx:
        e = ents[i]
        for li in range(len(e["loads"])):
            tasks.append((e, "cell", li, level, work, of, it, timeout_min * 60, generate_only, None))
        if grid_every and i % grid_every == 0 and level == "coarse":
            tasks.append((e, "cell", 0, "medium", work, of, it, timeout_min * 60, generate_only, None))
    rows = {i: dict(pos=i, idx=ents[i]["idx"], why=ents[i]["why"], cases={}) for i in idx}
    packs = os.path.join(out, f"cases_s{shard:03d}")
    os.makedirs(packs, exist_ok=True)

    def collect(pos, res, case):
        rows[pos]["cases"][res["name"]] = res
        if os.path.isdir(case):
            _pack(case, os.path.join(packs, f"v{pos:04d}_i{rows[pos]['idx']}_{res['name']}.tar.gz"))
            if not keep_work:
                shutil.rmtree(case, ignore_errors=True)
        print(f"[{shard}/{shards}] вариант {pos} {res['name']}: {res.get('run_status')} за {res['wall_s']} с", flush=True)
    liners = []
    with ProcessPoolExecutor(max_workers=max(1, workers)) as ex:
        futs = {ex.submit(task, t): t for t in tasks}
        while futs:
            for f in as_completed(list(futs)):
                t = futs.pop(f)
                pos, res, case = f.result()
                collect(pos, res, case)
                if (t[1] == "cell" and t[2] == 0 and t[3] == level and res.get("run_status") == "ok"
                        and res["result"].get("outlet", {}).get("m", 0) > 0):
                    e = ents[pos]
                    lt = (e, "liner", 0, level, work, of, it, timeout_min * 60, generate_only,
                          cell_outlet_state(res, e))
                    futs[ex.submit(task, lt)] = lt
                    liners.append(pos)
                elif t[1] == "cell" and t[2] == 0 and t[3] == level and generate_only:
                    e = ents[pos]
                    fake = dict(T=e["loads"][0]["hot_init"]["T_eq"], Y=e["loads"][0]["hot_init"]["Y"],
                                m_dot_can=(e["concept"]["a_pz"] * e["loads"][0]["W_air"] + e["loads"][0]["W_fuel"])
                                / e["n_can"])
                    lt = (e, "liner", 0, level, work, of, it, timeout_min * 60, True, fake)
                    futs[ex.submit(task, lt)] = lt
                break
    for i in idx:
        st = [c.get("run_status") for c in rows[i]["cases"].values()]
        rows[i]["status"] = "ok" if st and all(s in ("ok", "только сетка") for s in st) else \
            ("частично" if any(s == "ok" for s in st) else "ошибка")
    of_ver = os.environ.get("WM_PROJECT_VERSION", "")
    res = dict(schema=SCHEMA, kv=data["kv"], n=n, shard=shard, shards=shards, n_rows=len(rows),
               inputs_fingerprint=data["fingerprint"], inputs_file=os.path.basename(inputs), code_hash=code_hash(),
               level=level, iters=it, generate_only=generate_only,
               method=("уровень 1 воронки камеры: comb_cfd (осесимметричный клин, reactingFoam LTS, k-epsilon, EDC v2005, "
                       "GRI-Mech 3.0 + TDAC, P1) + comb_post (итоги без function objects)"),
               software=f"OpenFOAM {of_ver or 'не определено'}, python {sys.version.split()[0]}, numpy {CP.np.__version__}",
               date=CF.DATE, host=dict(cpu=os.cpu_count(), workers=workers, machine=platform.machine()),
               runtime_s=round(time.time() - t0, 1), rows=[rows[i] for i in idx])
    path = os.path.join(out, shard_name(data["kv"], n, shard, shards))
    with open(path + ".tmp", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, separators=(",", ":"), default=float)
    os.replace(path + ".tmp", path)
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    print(f"записано {path}: {len(rows)} вариантов за {res['runtime_s']} с", flush=True)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--n", type=int, required=True, help="сколько вариантов входного файла считать (первые n)")
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True, help="число частей (на AWS 10-50)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--level", default="coarse", choices=list(CF.MESH))
    ap.add_argument("--of", default="2512", help="версия словарей: 2512 (по умолчанию) или 1912 (местная обкатка)")
    ap.add_argument("--iters", default="", help="без химии, с химией, шаг записи: например 600,3000,500")
    ap.add_argument("--grid-every", type=int, default=10, help="каждый k-й вариант - ещё сетка medium (0 - нет)")
    ap.add_argument("--timeout-min", type=int, default=180, help="сторож одного кейса, мин")
    ap.add_argument("--generate-only", action="store_true")
    ap.add_argument("--keep-work", action="store_true")
    a = ap.parse_args(argv)
    it = None
    if a.iters:
        c, h, w = (int(x) for x in a.iters.split(","))
        it = dict(cold=c, hot=h, write=w)
    if shutil.which("reactingFoam") is None or shutil.which("blockMesh") is None:
        sys.exit("нет OpenFOAM в окружении (reactingFoam, blockMesh)")
    run(a.inputs, a.n, a.shard, a.shards, a.out, a.workers, a.level, a.of, it, a.grid_every, a.timeout_min,
        a.generate_only, a.keep_work)


if __name__ == "__main__":
    main()
