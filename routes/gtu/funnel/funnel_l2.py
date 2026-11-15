"""Уровень 2 воронки двигателя-заменителя: геометрия колёс, грубый CFD ступеней (SU2), осесимметричная прочность всех колёс и
дисков (CalculiX) для 50 разных уровня 1. Отбор на уровень 3 - 5-10 лучших прошедших.

  python3 funnel_l2.py probe          # короткая локальная проба -> l2_cases_probe_KV-3.json, funnel_l2_probe_results.json
  python3 funnel_l2.py prepare --input funnel_l1_top50.json [--top 50] --out l2_cases.json
  python3 funnel_l2.py run --cases l2_cases.json --what geom,fea --shards 10 --shard 0,5 --workers "$NPROC" --out "$GTU_OUT"
  python3 funnel_l2.py run --cases l2_cases.json --what su2 --level coarse --shards 17 --shard 3 --job-max-min 420 --out "$GTU_OUT"
  python3 funnel_l2.py merge --cases l2_cases.json --in <папка частей> --what geom,fea[,su2] --shards 10 --out-dir <папка>

Узлы платформы (не дублируются): геометрия и сетка канала центробежного колеса - ../platform/geom/impeller.py (--case custom);
линия частоты в SU2 - ../platform/su2/char_points.py run (база - стратегия разгона sweep_ramp.py или план smoke075, точки - рестарт
от базы); осесимметричный диск - ../platform/fea/disk.py через strength_l1.run_disk_fea. Деление на части - одно, по номеру
варианта в порядке отбора уровня 1 (div_rank): вариант r идёт в часть r % S, все его колёса и диски - в ту же часть.
Чего нет и что помечено как ожидающее (pending): CFD турбин и осевых/диагональных ступеней компрессора (нет генератора геометрии
решётки), одномерная камера сгорания (стык с воронкой камеры ../combustor/, поля funnel_interface). Всё - оценка.
AWS не вызывается; локально (4 ядра) - только prepare, geom с --no-mesh, fea и --dry для su2.
"""
import argparse, glob, hashlib, json, math, os, shlex, subprocess, sys, tempfile, time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
REL = "routes/gtu"
IMPELLER = os.path.join(GTU, "platform", "geom", "impeller.py")
CHAR = os.path.join(GTU, "platform", "su2", "char_points.py")
sys.path.insert(0, HERE)
import strength_l1 as ST          # noqa: E402

DATE = "2026-11-14"
T0_ISO, P0_ISO, GAMMA = 288.15, 101325.0, 1.4
MESH_L2 = dict(nr=80, nz=12)          # уровень 1 - 40×6 (strength_l1.disk_case); уровень 2 - вдвое мельче по каждой оси
T_HUB_DROP = 200.0                    # диск осевой турбины: ступица холоднее обода на 200 K по квадрату радиуса (как frame3_fea, оценка)
R_BORE = 0.20                         # радиус отверстия диска осевой ступени / радиус обода (оценка)
PFRACS = {"coarse": [0.6, 0.75, 0.9, 1.05], "medium": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]}
POINT_ITER = {"coarse": 3000, "medium": 6000}
# критерии прохода уровня 2 (README, «Уровни 2-3»; оценки)
PASS = dict(burst_120_min=1.0, life_h_min=25000.0, d_eta_su2_max=0.02, d_pi_su2_rel_max=0.05, mesh_ok=True)


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


# =========================================================================== prepare
def _stage_inlets(ml):
    """полные T и p на входе каждой ступени компрессора по степеням и политропным КПД одномерной модели (без потерь между
    ступенями - оценка)"""
    T, p, out = T0_ISO, P0_ISO, []
    for st in ml:
        out.append((T, p))
        pr, ep = float(st["PRs"]), float(st.get("eta_p") or 0.85)
        T = T * pr ** ((GAMMA - 1) / (GAMMA * ep))
        p = p * pr
    return out


def _axial_disc_case(name, disc, rpm, T_rim_C):
    """диск осевой ступени - равнопрочный (Стодола) по strength_l1.disc_axial: t(r) = t_об·exp(ρω²(R² − r²)/(2σ_d)), на ободе -
    растягивающая нагрузка лопаток p = σ_d (так выбрана толщина обода); материал - IN-713C (как уровень 1)"""
    R, s_d = float(disc["R"]), float(disc["sigma_d_MPa"]) * 1e6
    rho = float(ST._data()["mats"]["in713c"]["density"])
    om = 2 * math.pi * float(rpm) / 60
    rs = [R * (R_BORE + (1 - R_BORE) * k / 10) for k in range(11)]
    t = [float(disc["t_rim"]) * math.exp(rho * om ** 2 * (R ** 2 - r ** 2) / (2 * s_d)) for r in rs]
    r_in = rs[0]
    T = [T_rim_C - T_HUB_DROP * (1 - ((r - r_in) / (R - r_in)) ** 2) for r in rs]
    return dict(name=name, material="in713c", rpm=float(rpm), profile=dict(r=rs, t=t, hb=[0.0] * len(rs)), profile_interp="pchip",
                blade_fill=0.0, T_C=dict(r=rs, T=T), burst_k=1.0, rim=dict(p_MPa=s_d / 1e6), mesh=dict(MESH_L2),
                note="уровень 2: равнопрочный диск по strength_l1.disc_axial, обод - температура лопатки (верхняя граница), ступица -200 K (оценка)")


def variant_case(v):
    """один вариант отбора уровня 1 -> задания уровня 2"""
    key = f"r{int(v['div_rank']):02d}_j{v['j']}"
    ml, geo, strg, c = v["meanline"]["compressor"], v["geometry"]["compressor"], v["strength"], v["concept"]
    cores = int(c.get("cores") or 1)
    G_core = float(v["G"]) / cores
    inl = _stage_inlets(ml)
    stages, pending, disks = [], [], []
    solid = c["rotor_join"] in ("заодно", "сварка")
    for i, (st, g) in enumerate(zip(ml, geo)):
        sid = f"{key}_c{i + 1}"
        if g.get("type") != "центробежная":
            pending.append(dict(id=sid, what=f"CFD ступени «{g.get('type')}»", why="нет генератора геометрии этой ступени (impeller.py - "
                                                                                   "только центробежное колесо)"))
            continue
        T01, p01 = inl[i]
        args = dict(m=round(G_core, 4), pr=round(float(st["PRs"]), 5), rpm=round(float(v["n_rpm"]), 1), d2=round(float(st["D2"]), 5),
                    z=int(round(float(st.get("Zfull") or g.get("blades")))), b2=round(float(st["b2"]), 5), beta2b=round(float(st["beta2b_deg"]), 3),
                    d1s=round(float(st["D1t"]), 5), nu=round(float(st["D1h"]) / float(st["D1t"]), 4), gap=float(st.get("cl") or 3e-4),
                    t01=round(T01, 2), p01=round(p01, 1))
        if not g.get("splitters"):
            args["no_splitters"] = True
        if st.get("dbl"):
            args["double_sided"] = True
        stages.append(dict(id=sid, stage=i + 1, stage_n=len(stages) + 1, args=args, meanline=dict(PRs=st["PRs"], eta_p=st.get("eta_p"), M1t=st.get("M1t"),
                                                                          U2=st.get("U2"), phi=st.get("phi"))))
    for w, st in zip(strg.get("comp", []), ml):
        if w.get("material_data") in ("ti64", "in713c"):
            T_in = float(w.get("T_in_C") if w.get("T_in_C") is not None else 15.0)
            case = ST.disk_case(f"l2_{key}_c{w['stage']}", "comp", w["material_data"], v["n_rpm"], st["D2"], T_in, w["T_out_C"], solid)
            case["mesh"] = dict(MESH_L2)
            disks.append(dict(id=f"{key}:comp{w['stage']}", case=case))
        else:
            pending.append(dict(id=f"{key}:comp{w['stage']}", what="прочность колеса компрессора",
                                why=f"материал {w.get('material_data')} без данных в materials.json"))
    t = strg.get("turb") or {}
    if c["turb"] == "Р" and c.get("mat_tw") != "Si3N4" and t:
        D4 = v["meanline"]["turbine"].get("D4")
        if D4:
            case = ST.disk_case(f"l2_{key}_t", "rturb", "in713c", v["n_rpm"], D4, t["T_rim_C"] - 240.0, t["T_rim_C"], True)
            case["mesh"] = dict(MESH_L2)
            disks.append(dict(id=f"{key}:turb", case=case))
    elif c["turb"] == "О" and t.get("disc"):
        disks.append(dict(id=f"{key}:turb_disc", case=_axial_disc_case(f"l2_{key}_td", t["disc"], v["n_rpm"], float(t["T_blade_C"]))))
    if c.get("mat_tw") == "Si3N4":
        pending.append(dict(id=f"{key}:turb", what="прочность керамического колеса", why="disk.py без хрупкого критерия; Вейбулл - уровень 1, 3D - уровень 3"))
    pt = strg.get("pt") or {}
    if pt.get("disc") and v["geometry"].get("power_turbine"):
        disks.append(dict(id=f"{key}:pt_disc", case=_axial_disc_case(f"l2_{key}_ptd", pt["disc"], v["geometry"]["power_turbine"]["n_rpm"],
                                                                      float(pt["T_blade_C"]))))
    nt = (v["meanline"]["turbine"].get("stages") or [None])
    pending.append(dict(id=f"{key}:turbine_cfd", what=f"CFD турбины ({c['turb']}, {len(nt)} ст.) и силовой турбины",
                        why="нет генератора геометрии турбинной решётки и сетки канала (узел geom - только центробежное колесо)"))
    cb = v["geometry"].get("combustor") or {}
    comb = dict(type=cb.get("type"), D=cb.get("D"), L=cb.get("L"), injectors=cb.get("injectors"), pi=v.get("pi"), t41_C=v.get("t41_C"),
                G=v.get("G"), cores=cores, note="одномерная камера - стык с воронкой камеры (../combustor, merge_l1 -> feed_to_engine_funnel, "
                                                "targets.json funnel_interface): потери, NOx/CO, RTDF, температура стенки по типу камеры и p3/T3/T41")
    pending.append(dict(id=f"{key}:comb1d", what="одномерная камера сгорания", why="стык с воронкой камеры не написан (данные - её merge_l1)"))
    return dict(key=key, div_rank=v["div_rank"], j=v["j"], idx=v["idx"], arch=v["arch"], niche=v["niche"], scheme=v.get("scheme"),
                score_l1=v["score"], eta_th_shaft_pct=v.get("eta_th_shaft_pct"), n_rpm=v["n_rpm"], stages=stages, disks=disks, comb=comb,
                pending=pending)


def prepare(input_path, top, out, select=None):
    """select - funnel_l2_results.json: только варианты l3_select (уровень 3), в его порядке; номер части = место в списке (slot)"""
    d = json.load(open(input_path, encoding="utf-8"))
    sel = sorted(d["selected"], key=lambda v: v["div_rank"])
    if select:
        keys = json.load(open(select, encoding="utf-8"))["l3_select"]
        by = {f"r{int(v['div_rank']):02d}_j{v['j']}": v for v in sel}
        vs = [dict(variant_case(by[k]), slot=i) for i, k in enumerate(keys) if k in by][:top]
    else:
        vs = [variant_case(v) for v in sel[:top]]
    doc = dict(title="уровень 2: задания по вариантам отбора уровня 1", date=DATE, kv=d.get("kv"), status="оценка",
               method="колёса - impeller.py --case custom по размерам одномерной модели уровня 1; диски - disk.py (сетка "
                      f"{MESH_L2['nr']}×{MESH_L2['nz']}); вход ступени - по степеням и политропным КПД (оценка)",
               input=dict(file=os.path.basename(input_path), sha256=_sha(input_path), code_sha_l1=d.get("code_sha"), M=d.get("M")),
               shard_rule=("вариант на месте slot списка l3_select -> часть slot % S" if select else
                           "вариант с номером r в порядке div_rank -> часть (r - 1) % S (одно деление)"),
               select=(dict(file=os.path.basename(select), sha256=_sha(select)) if select else None),
               n_variants=len(vs), n_stages=sum(len(x["stages"]) for x in vs), n_disks=sum(len(x["disks"]) for x in vs),
               n_pending=Counter(p["what"].split(" (")[0] for x in vs for p in x["pending"]), variants=vs)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)
    return doc


# =========================================================================== run
def _impeller_cmd(st, level, out_dir, threads, no_mesh=False):
    a = st["args"]
    cmd = [sys.executable, IMPELLER, "--case", "custom", "--level", level, "--out", out_dir, "--threads", str(threads)]
    for k in ("m", "pr", "rpm", "d2", "z", "b2", "beta2b", "d1s", "nu", "gap", "t01", "p01"):
        cmd += [f"--{k}", str(a[k])]
    if a.get("no_splitters"):
        cmd.append("--no-splitters")
    if a.get("double_sided"):
        cmd.append("--double-sided")
    if no_mesh:
        cmd.append("--no-mesh")
    return cmd


def _impeller_summary(out_dir):
    js = sorted(glob.glob(os.path.join(out_dir, "impeller_custom*.json")))
    if not js:
        return None
    r = json.load(open(js[-1], encoding="utf-8"))
    d1, ar, mesh = r.get("d1") or {}, r.get("area") or {}, r.get("mesh") or {}
    return dict(U2=d1.get("U2"), b2=d1.get("b2"), M1s_rel=d1.get("M1s_rel"), beta2b=d1.get("beta2b"), slip=d1.get("slip_wiesner"),
                A_min_over_A1=ar.get("min_over_A1"), warnings=r.get("warnings"), mesh_nodes=(mesh.get("quality") or {}).get("n_nodes"),
                mesh_quality=({k: v for k, v in (mesh.get("quality") or {}).items() if k != "cell_types"} or None), t_total_s=r.get("t_total_s"))


def run_geom(stages, work, workers, level="coarse", no_mesh=False):
    def one(st):
        od = os.path.join(work, st["id"])
        os.makedirs(od, exist_ok=True)
        t = time.time()
        p = subprocess.run(_impeller_cmd(st, level, od, 1, no_mesh), capture_output=True, text=True)
        return st["id"], dict(rc=p.returncode, seconds=round(time.time() - t, 1), level=level, no_mesh=no_mesh,
                              summary=_impeller_summary(od), stderr_tail=(p.stderr[-600:] if p.returncode else ""))
    with ThreadPoolExecutor(max(1, workers)) as ex:
        return dict(ex.map(one, stages))


def run_fea(disks, work, workers):
    def one(it):
        try:
            return it["id"], ST.run_disk_fea(it["case"], os.path.join(work, it["case"]["name"]))
        except Exception as e:                           # noqa: BLE001 - отказ ccx записывается, часть не рушится
            return it["id"], dict(error=f"{type(e).__name__}: {e}")
    with ThreadPoolExecutor(max(1, workers)) as ex:
        return dict(ex.map(one, disks))


def su2_commands(st, level, work, out, nproc, job_max_min, base_strategy=None, n_stage_left=1, pfracs=None):
    """команды одной ступени: сетка impeller.py (все ядра) -> char_points.py run (база + точки линии частоты на сетке level)"""
    geo = os.path.join(work, st["id"] + "_" + level)
    pf = pfracs or PFRACS[level]
    npp = max(1, int(nproc) // len(pf))
    budget = max(20, int(job_max_min / max(1, n_stage_left)))
    c1 = _impeller_cmd(st, level, geo, nproc)
    c2 = [sys.executable, CHAR, "run", "--case", "custom", "--level", level, "--geom-out", geo, "--label", st["id"],
          "--pfracs", ",".join(map(str, pf)), "--base-pfrac", "0.75", "--np-per-point", str(npp), "--point-iter", str(POINT_ITER[level]),
          "--point-pramp", "300", "--time-limit", str(budget * 60), "--work", os.path.join(work, "runs_" + st["id"]),
          "--out", os.path.join(out, st["id"])]
    c2 += (["--base-strategy", base_strategy] if base_strategy else ["--base-plan", "smoke075"])
    return [c1, c2]


def run_su2(stages, work, out, nproc, job_max_min, level, base_strategy, dry, pfracs=None):
    res, t_start = {}, time.time()
    for k, st in enumerate(stages):
        left = max(1.0, job_max_min - (time.time() - t_start) / 60 - 10)
        cmds = su2_commands(st, level, work, out, nproc, left, base_strategy, len(stages) - k, pfracs)
        if dry:
            res[st["id"]] = dict(dry=True, commands=[" ".join(shlex.quote(x) for x in c) for c in cmds])
            continue
        rc, t = [], time.time()
        for c in cmds:
            p = subprocess.run(c, capture_output=True, text=True)
            rc.append(p.returncode)
            if p.returncode:
                break
        summ = sorted(glob.glob(os.path.join(out, st["id"], "**", "char_summary.json"), recursive=True))
        res[st["id"]] = dict(rc=rc, seconds=round(time.time() - t, 1), level=level, char_summary=(os.path.relpath(summ[-1], out) if summ else None))
    return res


def mine(doc, shards, shard_list):
    return [v for v in doc["variants"] if v.get("slot", v["div_rank"] - 1) % shards in shard_list]


def parse_shards(spec, shards):
    if spec in ("all", "*"):
        return list(range(shards))
    out = set()
    for part in str(spec).split(","):
        if "-" in part:
            a, b = part.split("-"); out.update(range(int(a), int(b) + 1))
        elif part.strip():
            out.add(int(part))
    bad = [k for k in out if not 0 <= k < shards]
    if bad:
        raise ValueError(f"части {bad} вне 0..{shards - 1}")
    return sorted(out)


def run(a):
    doc = json.load(open(a.cases, encoding="utf-8"))
    whats = a.what.split(",")
    if "su2" in whats and a.base_strategy and not a.dry and not os.path.exists(a.base_strategy):
        # заглушка plan_sweep_WINNER.json до вердикта перебора разгона (ворота jobs_l2.json): отказ сразу, а не после сетки и базы
        raise SystemExit(f"нет стратегии разгона {a.base_strategy}: подставить план-победитель перебора (make_jobs_l23.py --sweep-plan)")
    work = a.work or tempfile.mkdtemp(prefix="l2_")
    os.makedirs(a.out, exist_ok=True)
    for k in parse_shards(a.shard, a.shards):
        vs = mine(doc, a.shards, [k])
        part = dict(level=2, cases=os.path.basename(a.cases), cases_sha256=_sha(a.cases), shard=k, shards=a.shards,
                    variants=[v["key"] for v in vs], date=DATE, nproc=a.workers, results={})
        t = time.time()
        if "geom" in whats:
            part["results"]["geom"] = run_geom([s for v in vs for s in v["stages"]], work, a.workers, a.level, a.no_mesh)
        if "fea" in whats:
            part["results"]["fea"] = run_fea([d for v in vs for d in v["disks"]], work, a.workers)
        if "su2" in whats:
            sts = [s for v in vs for s in v["stages"] if a.stage is None or s["stage_n"] == a.stage]
            pf = [float(x) for x in a.pfracs.split(",")] if a.pfracs else None
            part["results"]["su2"] = run_su2(sts, work, a.out, a.workers, a.job_max_min, a.level, a.base_strategy, a.dry, pf)
            part.update(stage=a.stage, pfracs=pf, tag=a.tag)
        part["seconds"] = round(time.time() - t, 1)
        tag = (f"_st{a.stage}" if a.stage else "") + (f"_{a.tag}" if a.tag else "")
        name = f"funnel_l2_{'-'.join(whats)}_{a.level}{tag}_s{k:02d}of{a.shards}.json"
        with open(os.path.join(a.out, name), "w", encoding="utf-8") as fh:
            json.dump(part, fh, ensure_ascii=False, indent=1)
        print(os.path.join(a.out, name), part["seconds"], "с")
    return 0


# =========================================================================== merge
def merge(a):
    doc = json.load(open(a.cases, encoding="utf-8"))
    whats = a.what.split(",")
    parts = {}            # ключ части: (что, номер части, ступень, метка половины точек) - у su2 уровня 3 на часть несколько задач
    for p in sorted(glob.glob(os.path.join(a.inp, "**", "funnel_l2_*.json"), recursive=True)):
        d = json.load(open(p, encoding="utf-8"))
        if d.get("level") != 2 or d.get("shards") != a.shards:
            continue
        for w in whats:
            if w in d["results"]:
                pk = (w, d["shard"], d.get("stage") if w == "su2" else None, d.get("tag", "") if w == "su2" else "")
                if pk in parts and parts[pk]["results"][w] != d["results"][w]:
                    raise SystemExit(f"часть {d['shard']} ({w}) задана дважды с разным содержанием: {p}")
                if d["cases_sha256"] != _sha(a.cases):
                    raise SystemExit(f"часть {p} посчитана по другому файлу заданий")
                parts[pk] = d
    problems = []
    for w in whats:
        have = {pk[1] for pk in parts if pk[0] == w}
        miss = [k for k in range(a.shards) if k not in have]
        if miss:
            problems.append(f"{w}: нет частей {miss}")
    got = {w: {} for w in whats}
    for pk, d in parts.items():
        w = pk[0]
        for key, r in d["results"][w].items():
            if w == "su2":                    # половины точек одной ступени (метки a, b) складываются
                got[w].setdefault(key, []).append(dict(r, tag=pk[3]))
                continue
            if key in got[w]:
                problems.append(f"{w}: {key} задвоен")
            got[w][key] = r
    rows = []
    for v in doc["variants"]:
        exp = dict(geom=[s["id"] for s in v["stages"]], su2=[s["id"] for s in v["stages"]], fea=[d["id"] for d in v["disks"]])
        r = dict(key=v["key"], arch=v["arch"], niche=v["niche"], score_l1=v["score_l1"], eta_th_shaft_pct=v["eta_th_shaft_pct"], fails=[])
        for w in whats:
            lost = [x for x in exp[w] if x not in got[w]]
            if lost:
                problems.append(f"{w}: потеряны {lost}")
        if "geom" in whats:
            g = [got["geom"].get(x) or {} for x in exp["geom"]]
            r["geom_ok"] = all(x.get("rc") == 0 for x in g)
            if not r["geom_ok"]:
                r["fails"].append("геометрия колеса не построена")
        if "fea" in whats:
            f = [got["fea"].get(x) or {} for x in exp["fea"]]
            bm = [x.get("burst_margin_120") for x in f if x.get("burst_margin_120") is not None]
            life = [x.get("life_local_h") for x in f if x.get("life_local_h") is not None]
            r["fea"] = {x.split(":", 1)[1]: {k: (round(y[k], 3) if isinstance(y.get(k), float) and y[k] < 1e9 else y.get(k))
                                              for k in ("burst_margin_120", "life_local_h", "vm_max_MPa", "error") if k in y}
                        for x, y in zip(exp["fea"], f)}
            r["burst_120_min"] = min(bm) if bm else None
            r["life_h_min"] = min(life) if life else None
            if any("error" in x for x in f) or len(bm) < len(f):
                r["fails"].append("CalculiX не посчитал часть дисков")
            if bm and min(bm) < PASS["burst_120_min"]:
                r["fails"].append(f"запас по разрушению при 120 % {min(bm):.2f} < {PASS['burst_120_min']}")
            if life and min(life) < PASS["life_h_min"]:
                r["fails"].append(f"ресурс по ползучести {min(life):.0f} ч < {PASS['life_h_min']:.0f}")
        if "su2" in whats:
            s = [got["su2"].get(x) or [] for x in exp["su2"]]
            r["su2_done"] = sum(1 for x in s if any(y.get("char_summary") for y in x))
            if r["su2_done"] < len(s):
                r["fails"].append("SU2: линия частоты не получена у части ступеней (сравнение с одномерной моделью - вручную по char_summary)")
        r["pass"] = not r["fails"]
        rows.append(r)
    ok = sorted([r for r in rows if r["pass"]], key=lambda r: -r["score_l1"])
    l3, per_arch = [], Counter()
    for r in ok:
        if per_arch[r["arch"]] < 2 and len(l3) < 10:
            l3.append(r["key"]); per_arch[r["arch"]] += 1
    os.makedirs(a.out_dir, exist_ok=True)
    no_mesh = any(r.get("no_mesh") for pk, d in parts.items() if pk[0] == "geom" for r in d["results"]["geom"].values())
    out = dict(title="уровень 2: итог", date=DATE, status="оценка" if not problems else "неполно",
               method=("колёса - impeller.py --case custom (STEP" + (", без сетки - локальная проба" if no_mesh else " и сетка канала") +
                       "); диски - disk.py (CalculiX, осесимметричная модель, сетка " + f"{MESH_L2['nr']}×{MESH_L2['nz']}" +
                       "; диски осевых ступеней - равнопрочный профиль, обод при температуре лопатки - верхняя граница); SU2 - "
                       "char_points.py (если задано). Всё - оценка"),
               software=f"python {sys.version.split()[0]}; impeller.py, disk.py (ccx), char_points.py узлов платформы",
               inputs=[dict(file=os.path.basename(a.cases), sha256=_sha(a.cases), source=doc.get("input"))]
               + [dict(part=f"{pk[0]}:s{pk[1]}" + (f":st{pk[2]}" if pk[2] else "") + (f":{pk[3]}" if pk[3] else ""), seconds=d.get("seconds"),
                       nproc=d.get("nproc")) for pk, d in sorted(parts.items(), key=lambda x: str(x[0]))],
               whats=whats, shards=a.shards,
               completeness=dict(ok=not problems, problems=problems), pass_rule=PASS, n=len(rows), n_pass=len(ok),
               l3_select=l3, l3_rule="прошедшие по убыванию оценки уровня 1, не больше 2 на архетип, до 10 (не меньше 5 - иначе точка остановки)",
               rows=rows)
    path = os.path.join(a.out_dir, "funnel_l2_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print(path, "полнота:", "да" if not problems else problems[:5], "прошли:", len(ok), "на уровень 3:", len(l3))
    return 0 if not problems else 2


def probe(workers=4):
    """короткая локальная проба (решение владельца: локально - только пробы): отобранные пробой уровня 1 (funnel_l1_top50.json, до 10)
    -> l2_cases_probe_KV-3.json -> геометрия колёс без сетки (STEP) и disk.py, 2 части -> слияние -> funnel_l2_probe_results.json"""
    cases = os.path.join(HERE, "l2_cases_probe_KV-3.json")
    prepare(os.path.join(HERE, "funnel_l1_top50.json"), 10, cases)
    tmp = tempfile.mkdtemp(prefix="l2_probe_")
    try:
        run(argparse.Namespace(cases=cases, what="geom,fea", shards=2, shard="all", workers=workers, level="coarse", no_mesh=True,
                               base_strategy=None, job_max_min=60.0, dry=False, stage=None, pfracs=None, tag="",
                               work=os.path.join(tmp, "w"), out=os.path.join(tmp, "parts")))
        rc = merge(argparse.Namespace(cases=cases, inp=os.path.join(tmp, "parts"), what="geom,fea", shards=2, out_dir=tmp))
        r = json.load(open(os.path.join(tmp, "funnel_l2_results.json"), encoding="utf-8"))
        r["title"] = "уровень 2: короткая локальная проба (отобранные пробой уровня 1; геометрия без сетки, disk.py)"
        r["probe_note"] = ("локальная проба на 4 ядрах: колёса - только STEP (--no-mesh), SU2 не считался; полный уровень 2 - задачи AWS "
                           "jobs_l2.json по 50 разным слитого уровня 1")
        with open(os.path.join(HERE, "funnel_l2_probe_results.json"), "w", encoding="utf-8") as fh:
            json.dump(r, fh, ensure_ascii=False, indent=1)
        return rc
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("probe")
    p = sp.add_parser("prepare")
    p.add_argument("--input", required=True)
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--out", required=True)
    p.add_argument("--select", default=None, help="funnel_l2_results.json: задания уровня 3 - только l3_select")
    r = sp.add_parser("run")
    r.add_argument("--cases", required=True)
    r.add_argument("--what", default="geom,fea")
    r.add_argument("--shards", type=int, required=True)
    r.add_argument("--shard", default="all")
    r.add_argument("--workers", type=int, default=int(os.environ.get("NPROC") or os.cpu_count() or 1))
    r.add_argument("--level", default="coarse", choices=list(PFRACS))
    r.add_argument("--no-mesh", action="store_true", help="geom: только STEP (локальная проба)")
    r.add_argument("--base-strategy", default=None, help="su2: стратегия разгона plan_sweep_*.json (после вердикта перебора)")
    r.add_argument("--job-max-min", type=float, default=600.0)
    r.add_argument("--dry", action="store_true", help="su2: только команды")
    r.add_argument("--stage", type=int, default=None, help="su2: только центробежная ступень с этим номером (1, 2) у каждого варианта")
    r.add_argument("--pfracs", default=None, help="su2: доли давления на выходе через запятую (иначе PFRACS уровня сетки)")
    r.add_argument("--tag", default="", help="метка в имени файла части (половины точек одной ступени: a, b)")
    r.add_argument("--work", default=None)
    r.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    m = sp.add_parser("merge")
    m.add_argument("--cases", required=True)
    m.add_argument("--in", dest="inp", required=True)
    m.add_argument("--what", default="geom,fea")
    m.add_argument("--shards", type=int, required=True)
    m.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        d = prepare(a.input, a.top, a.out, a.select)
        print(a.out, d["n_variants"], "вариантов,", d["n_stages"], "колёс,", d["n_disks"], "дисков; ожидают:", dict(d["n_pending"]))
        return 0
    if a.cmd == "run":
        return run(a)
    if a.cmd == "probe":
        return probe()
    return merge(a)


if __name__ == "__main__":
    sys.exit(main())
