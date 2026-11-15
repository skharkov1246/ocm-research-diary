"""Воронка камеры сгорания KV, уровни 2 и 3: одна часть - трёхмерные секторы элемента (OpenFOAM, параллельно).

Уровень 2 (--level L2), на вариант (comb_cfd3d.build_element3d_case, RANS LTS, k-epsilon, EDC, GRI-Mech 3.0 + TDAC):
  * элемент 3D на 100 / 50 / 30 %;
  * запас по срыву: рестарты 100 % с топливом ×0,9; 0,8; 0,7; 0,6; 0,5 (расход топлива - скорость струи на входе
    последней записи); погасание - полнота < 0,9 или T_max < T3 + 300 K; запас = 1 / (последний горящий множитель);
  * запас по проскоку: рестарты с долей H2 в топливе 0,3; 0,5; 0,7; 1,0 по объёму при той же теплоте топлива;
    проскок - T в канале > T3 + 200 K; итог - наибольшая доля H2 без проскока.
Теплообмен стенки (сектор трубы chtMultiRegionFoam) и распыл дизеля (sprayFoam) - генератор сетки и словарей готов
частично (comb_cfd3d.py): в README и jobs_l2.json - «не готово к запуску» с перечнем того, что дописать.
Уровень 3 (--level L3): элемент 3D мельче, LES (WALE) нестационарно, датчики давления и интеграл тепловыделения;
с --forcing - широкополосное возбуждение расхода воздуха на входе для передаточной функции пламени
(comb_acoustic.identify_ftf).

Кейс - --np ядер (decomposePar, mpirun); одновременно --workers / --np кейсов. Итоги - comb_post.analyse.
Запуск:
  python comb_l23.py --level L2 --inputs l2_inputs_KV-3.json --n 40 --shard 3 --shards 20 --out "$GTU_OUT" --np 8
  python comb_l23.py --level L2 ... --generate-only     # только сетки (blockMesh, checkMesh)
"""
import argparse, json, math, os, shutil, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_cfd as CF    # noqa: E402
import comb_cfd3d as C3  # noqa: E402
import comb_post as CP   # noqa: E402
import comb_l1 as L1     # noqa: E402

SCHEMA = "kv.combustor.l23.shard/1"
LBO_STEPS = [0.9, 0.8, 0.7, 0.6, 0.5]
FB_STEPS = [0.3, 0.5, 0.7, 1.0]
ITERS2 = dict(cold=800, hot=3500, write=500, restart=1500)
LES = dict(flow_through=8, stats=6, max_co=0.5)
LES_LENGTH_FRAC = 0.35     # LES - канал и зона пламени: 35 % длины ячейки (смета, README «Уровень 3»)


def code_hash():
    import hashlib
    h = hashlib.sha256()
    for f in ("comb_cfd.py", "comb_cfd3d.py", "comb_post.py", "comb_l23.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def _decompose(case, n):
    CF.foam_file(os.path.join(case, "system", "decomposeParDict"), "dictionary", "decomposeParDict",
                 f"numberOfSubdomains {n};\n\nmethod          scotch;")


def _solve(case, log, n, timeout_s):
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    t0 = time.time()
    with open(os.path.join(case, log), "w") as f:
        try:
            if n > 1:
                subprocess.run(["decomposePar", "-force"], cwd=case, stdout=f, stderr=subprocess.STDOUT, timeout=600)
                r = subprocess.run(["mpirun", "-np", str(n), "--oversubscribe", "reactingFoam", "-parallel"], cwd=case,
                                   stdout=f, stderr=subprocess.STDOUT, timeout=timeout_s, env=env)
                subprocess.run(["reconstructPar", "-latestTime"], cwd=case, stdout=f, stderr=subprocess.STDOUT, timeout=600)
                for d in os.listdir(case):
                    if d.startswith("processor"):
                        shutil.rmtree(os.path.join(case, d), ignore_errors=True)
            else:
                r = subprocess.run(["reactingFoam"], cwd=case, stdout=f, stderr=subprocess.STDOUT, timeout=timeout_s)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = 124
    return rc, round(time.time() - t0, 1)


def _fd(case, path, entry, value):
    subprocess.run(["foamDictionary", path, "-entry", entry, "-set", value], cwd=case, capture_output=True)


def _end(case, n_add):
    t = CP.latest_times(case, 1)
    t0 = float(t[0]) if t else 0.0
    _fd(case, "system/controlDict", "endTime", str(int(t0 + n_add)))
    return t[0] if t else "0"


def _restart_fuel(case, scale=None, h2=None, desc=None):
    """следующий шаг цепочки: множитель расхода топлива или доля H2 по объёму при той же теплоте топлива"""
    t = CP.latest_times(case, 1)[0]
    j = desc["jet"]
    if scale is not None:
        u = j["u"] * scale
        Yc, Yh = None, None
    else:
        mwf = h2 * 2.016 + (1 - h2) * 16.043
        Yh = h2 * 2.016 / mwf
        Yc = 1 - Yh
        lhv = Yc * 50.0254e6 + Yh * 119.953e6                 # comb_bc_results.json → LHV_J_kg (CH4, H2)
        m = j["m"] * 50.0254e6 / lhv
        rho = desc["p3"] * mwf * 1e-3 / (CF.R_U * 288.15)
        u = j["u"] * (m / j["m"]) * (j["rho"] / rho)
    U0 = CP.read_field(os.path.join(case, t, "U"))["boundary"]["inlet_fuel"]["value"]
    k = u / max(math.sqrt(sum(x * x for x in U0)), 1e-30)
    _fd(case, f"{t}/U", "boundaryField.inlet_fuel.value", f"uniform ({U0[0] * k:.9g} {U0[1] * k:.9g} {U0[2] * k:.9g})")
    if Yc is not None:
        _fd(case, f"{t}/CH4", "boundaryField.inlet_fuel.value", f"uniform {Yc:.9g}")
        _fd(case, f"{t}/H2", "boundaryField.inlet_fuel.value", f"uniform {Yh:.9g}")


def _status_burn(a, T3):
    o = a.get("outlet", {})
    return (o.get("eta_c") or 0) >= 0.9 and (a.get("T_max") or 0) >= T3 + 300.0 and abs(a.get("mass_imbalance", 1)) < 0.05


def to_les(case, desc, forcing=False):
    """кейс элемента 3D -> LES (уровень 3): WALE, Эйлер, шаг по числу Куранта, датчики, интеграл тепловыделения"""
    CF.foam_file(os.path.join(case, "constant", "turbulenceProperties"), "dictionary", "turbulenceProperties",
                 "simulationType LES;\n\nLES\n{\n    LESModel        WALE;\n    turbulence      on;\n    printCoeffs     on;\n"
                 "    delta           cubeRootVol;\n}")
    s = open(os.path.join(case, "system", "fvSchemes")).read().replace("localEuler", "Euler")
    open(os.path.join(case, "system", "fvSchemes"), "w").write(s)
    zd, ze, rc = desc["z_dump"], desc["z_end"], desc["geometry"]["r_ch"]
    u_dump = desc["u_ch"] * (desc["geometry"]["r_ch"] / desc["geometry"]["R_c"]) ** 2   # средняя скорость за срезом
    t_ft = (ze - zd) / max(u_dump, 1.0)
    T_end = (LES["flow_through"] + LES["stats"]) * t_ft
    pts = " ".join(f"(0 0 {z:.9g})" for z in (zd - 0.5 * rc, zd + 2 * rc, zd + 6 * rc, zd + 0.5 * (ze - zd), ze - 1e-4))
    CF.foam_file(os.path.join(case, "system", "controlDict"), "dictionary", "controlDict",
                 f"application     reactingFoam;\nstartFrom       latestTime;\nstartTime       0;\nstopAt          endTime;\n"
                 f"endTime         {T_end:.9g};\ndeltaT          1e-8;\nadjustTimeStep  yes;\nmaxCo           {LES['max_co']};\n"
                 f"maxDeltaT       1e-5;\nwriteControl    adjustable;\nwriteInterval   {T_end / 4:.9g};\npurgeWrite      2;\n"
                 f"writeFormat     ascii;\nwritePrecision  9;\nwriteCompression off;\ntimeFormat      general;\ntimePrecision   10;\n"
                 f"runTimeModifiable true;\n\nfunctions\n{{\n    probes\n    {{\n        type            probes;\n"
                 f"        libs            (sampling);\n        fields          (p U T);\n        probeLocations  ({pts});\n    }}\n"
                 f"    heatRelease\n    {{\n        type            volFieldValue;\n        libs            (fieldFunctionObjects);\n"
                 f"        fields          (Qdot);\n        operation       volIntegrate;\n        regionType      all;\n"
                 f"        writeFields     false;\n    }}\n    average\n    {{\n        type            fieldAverage;\n"
                 f"        libs            (fieldFunctionObjects);\n        timeStart       {LES['flow_through'] * t_ft:.9g};\n"
                 f"        fields          (U {{ mean on; prime2Mean on; base time; }} T {{ mean on; prime2Mean on; base time; }}"
                 f" p {{ mean on; prime2Mean on; base time; }});\n    }}\n}}\n")
    fv = open(os.path.join(case, "system", "fvSolution")).read()
    fv = fv.replace("nCorrectors         2", "nCorrectors         3")
    open(os.path.join(case, "system", "fvSolution"), "w").write(fv)
    if forcing:
        U = CP.read_field(os.path.join(case, "0", "U"))["boundary"]["inlet_air"]["value"]
        n, rows = 4000, []
        fr = [50 + 25 * k for k in range(59)]                   # 50 ... 1500 Гц
        ph = [((k * 7919) % 360) * math.pi / 180 for k in range(len(fr))]   # фазы - детерминированно (Шрёдер не нужен)
        amp = 0.05 / math.sqrt(len(fr) / 2)
        for i in range(n + 1):
            t = T_end * i / n
            f = 1 + amp * sum(math.sin(2 * math.pi * fk * t + p) for fk, p in zip(fr, ph))
            rows.append(f"({t:.9g} ({U[0] * f:.9g} {U[1] * f:.9g} {U[2] * f:.9g}))")
        _fd(case, "0/U", "boundaryField.inlet_air.type", "uniformFixedValue")
        open(os.path.join(case, "constant", "forcingTable"), "w").write("(\n" + "\n".join(rows) + "\n)\n")
        _fd(case, "0/U", "boundaryField.inlet_air.uniformValue", "tableFile")
        _fd(case, "0/U", "boundaryField.inlet_air.uniformValueCoeffs", "{ file \"<constant>/forcingTable\"; }")
    return dict(T_end=T_end, t_flow_through=t_ft, forcing=forcing)


def task(args):
    e, li, level, work, of, np_, timeout_s, gen_only, kind = args
    t0 = time.time()
    name = f"{kind}_{int(round(100 * e['loads'][li]['load']))}_{level}"
    case = os.path.join(work, f"v{e['_pos']:04d}_i{e['idx']}", name)
    res = dict(name=name, kind=kind, load=e["loads"][li]["load"], level=level, steps=[])
    T3 = e["loads"][li]["T3_K"]
    inp = dict(T3=T3, lhv=e["lhv"])
    try:
        desc = C3.build_element3d_case(e, li, case, "L3" if level == "L3" else "L2", of, ITERS2,
                                       length_frac=LES_LENGTH_FRAC if level == "L3" else 1.0)
        res["n_cells"] = desc["n_cells"]
        for cmd in (["blockMesh"], ["checkMesh"]):
            rc, dt = L1._run(cmd, case, f"log.{cmd[0]}", 1800)
            res["steps"].append(dict(cmd=cmd[0], rc=rc, s=round(dt, 1)))
        if "Mesh OK." not in open(os.path.join(case, "log.checkMesh"), errors="replace").read():
            res["run_status"] = "checkMesh: сетка с ошибками"
            return e["_pos"], res, case
        if gen_only:
            res["run_status"] = "только сетка"
            res["wall_s"] = round(time.time() - t0, 1)
            return e["_pos"], res, case
        _decompose(case, np_)
        _fd(case, "constant/chemistryProperties", "chemistry", "off")
        _end(case, ITERS2["cold"])
        rc, dt = _solve(case, "log.reactingFoam.cold", np_, timeout_s * 0.3)
        res["steps"].append(dict(cmd="cold", rc=rc, s=dt))
        L1._run(["setFields"], case, "log.setFields", 600)
        _fd(case, "constant/chemistryProperties", "chemistry", "on")
        if level == "L3":
            res["les"] = to_les(case, desc, forcing=(kind == "les_forced"))
        else:
            _end(case, ITERS2["hot"])
        rc, dt = _solve(case, "log.reactingFoam.hot", np_, timeout_s * 0.5)
        res["steps"].append(dict(cmd="hot", rc=rc, s=dt))
        if rc != 0:
            res["run_status"] = "не хватило времени" if rc == 124 else "решатель упал"
            return e["_pos"], res, case
        res["result"] = CP.analyse(case, inp)
        res["run_status"] = "ok"
        if kind == "lbo":
            last_burning = 1.0
            for sc in LBO_STEPS:
                _restart_fuel(case, scale=sc, desc=desc)
                _end(case, ITERS2["restart"])
                rc, dt = _solve(case, f"log.lbo_{sc}", np_, timeout_s * 0.1)
                a = CP.analyse(case, inp)
                res["steps"].append(dict(cmd=f"lbo x{sc}", rc=rc, s=dt, burning=_status_burn(a, T3)))
                if rc != 0 or not _status_burn(a, T3):
                    break
                last_burning = sc
                desc["jet"]["u"] *= 1.0          # скорость следующего шага считается от исходной струи
            res["lbo_margin"] = 1.0 / last_burning
        if kind == "flashback":
            ok_h2 = 0.0
            for h2 in FB_STEPS:
                _restart_fuel(case, h2=h2, desc=desc)
                _end(case, ITERS2["restart"])
                rc, dt = _solve(case, f"log.fb_{h2}", np_, timeout_s * 0.1)
                a = CP.analyse(case, inp)
                fb = bool(a.get("flashback_or_autoignition"))
                res["steps"].append(dict(cmd=f"H2 {h2}", rc=rc, s=dt, flashback=fb))
                if rc != 0 or fb:
                    break
                ok_h2 = h2
            res["H2_vol_max"] = ok_h2
    except Exception as ex:
        res.update(run_status="ошибка", error=f"{type(ex).__name__}: {str(ex)[:300]}")
    res["wall_s"] = round(time.time() - t0, 1)
    return e["_pos"], res, case


def run(inputs, n, shard, shards, out, level, np_, workers, of="2512", timeout_min=600, generate_only=False):
    t0 = time.time()
    data = json.load(open(inputs, encoding="utf-8"))
    ents = data["entries"]
    work = os.path.join(out, f"work_{level}_s{shard:03d}")
    os.makedirs(work, exist_ok=True)
    tasks = []
    if level == "L2":                     # разбиение по варианту: i % shards == shard
        idx = L1.shard_indices(n, shard, shards)
        for i in idx:
            ents[i]["_pos"] = i
            tasks += [(ents[i], li, level, work, of, np_, timeout_min * 60, generate_only, "base") for li in (1, 2)]
            tasks += [(ents[i], 0, level, work, of, np_, timeout_min * 60, generate_only, k) for k in ("lbo", "flashback")]
    else:                                 # L3: разбиение по задаче (вариант × вид): t % shards == shard, t = 2·i + вид
        allt = [(i, k) for i in range(n) for k in ("les", "les_forced")]
        mine = [allt[t] for t in L1.shard_indices(len(allt), shard, shards)]
        idx = sorted(set(i for i, _ in mine))
        for i, k in mine:
            ents[i]["_pos"] = i
            tasks.append((ents[i], 0, level, work, of, np_, timeout_min * 60, generate_only, k))
    rows = {i: dict(pos=i, idx=ents[i]["idx"], cases={}) for i in idx}
    packs = os.path.join(out, f"cases_{level}_s{shard:03d}")
    os.makedirs(packs, exist_ok=True)
    with ProcessPoolExecutor(max_workers=max(1, workers // max(1, np_))) as ex:
        for f in as_completed([ex.submit(task, t) for t in tasks]):
            pos, res, case = f.result()
            rows[pos]["cases"][res["name"]] = res
            if os.path.isdir(case):
                L1._pack(case, os.path.join(packs, f"v{pos:04d}_{res['name']}.tar.gz"))
                shutil.rmtree(case, ignore_errors=True)
            print(f"[{level} {shard}/{shards}] {pos} {res['name']}: {res.get('run_status')} за {res.get('wall_s')} с", flush=True)
    res = dict(schema=SCHEMA, level=level, kv=data["kv"], n=n, shard=shard, shards=shards,
               inputs_fingerprint=data["fingerprint"], code_hash=code_hash(), np=np_, iters=ITERS2, les=LES,
               generate_only=generate_only, software=f"OpenFOAM {os.environ.get('WM_PROJECT_VERSION', 'не определено')}",
               date=CF.DATE, runtime_s=round(time.time() - t0, 1), rows=[rows[i] for i in idx])
    path = os.path.join(out, f"comb_{level.lower()}_{data['kv']}_n{n}_s{shard:03d}of{shards:03d}.json")
    json.dump(res, open(path, "w", encoding="utf-8"), ensure_ascii=False, default=float)
    shutil.rmtree(work, ignore_errors=True)
    print(f"записано {path}", flush=True)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--level", required=True, choices=["L2", "L3"])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--np", type=int, default=8, help="ядер на кейс")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--of", default="2512")
    ap.add_argument("--timeout-min", type=int, default=600)
    ap.add_argument("--generate-only", action="store_true")
    a = ap.parse_args(argv)
    if shutil.which("reactingFoam") is None:
        sys.exit("нет OpenFOAM в окружении")
    run(a.inputs, a.n, a.shard, a.shards, a.out, a.level, a.np, a.workers, a.of, a.timeout_min, a.generate_only)


if __name__ == "__main__":
    main()
