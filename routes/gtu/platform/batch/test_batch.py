"""Тесты узла batch: шаблоны AWS Batch, точка входа образа, смета, тарифы. Без сети и без AWS.

    python3 -m pytest -q routes/gtu/platform/batch/test_batch.py
"""
from __future__ import annotations

import copy
import csv
import io
import json
import os
import pathlib
import signal
import subprocess
import sys
import time
import tomllib

import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import cost_estimate  # noqa: E402
import entrypoint  # noqa: E402
import fetch_rates  # noqa: E402
import render  # noqa: E402

GIT = {"sha": "0" * 40, "branch": "test", "dirty": False}


def cfg():
    return render.load_config(HERE / "batch.toml")


def rates():
    return render.load_rates()


# ================================================================ render: шаблоны и сверка с моделью API
@pytest.mark.parametrize("region", ["us-east-1", "eu-north-1"])
@pytest.mark.parametrize("mode", ["line", "point"])
def test_render_passes_api_model_and_limits(region, mode):
    req, files, errs = render.render(cfg(), region, mode, rates(), GIT)
    assert errs == []
    sj = req["submit_array_job"]
    assert sj["arrayProperties"]["size"] == (10 if mode == "line" else 80)
    assert req["job_definition"]["containerProperties"]["image"].startswith(f"000000000000.dkr.ecr.{region}.")
    assert f":{region}:" in req["job_queue"]["computeEnvironmentOrder"][0]["computeEnvironment"]
    for text in files.values():
        assert "${" not in text.replace("${GTU_PREV_RESTART", "")   # bash-подстановка в команде — не шаблон


def test_render_key_fields():
    req, files, _ = render.render(cfg(), "us-east-1", "line", rates(), GIT)
    cr = req["compute_environment"]["computeResources"]
    assert cr["type"] == "SPOT" and cr["allocationStrategy"] == "SPOT_PRICE_CAPACITY_OPTIMIZED"
    assert cr["minvCpus"] == 0 and cr["maxvCpus"] + 64 <= 384
    assert cr["instanceTypes"] == ["c7a.16xlarge", "c7i.16xlarge", "c8i.16xlarge"]
    assert cr["tags"]["RunId"] == "gtu-EXAMPLE-kv3-char" and "Project" in cr["tags"]
    assert "bidPercentage" not in cr                     # пусто = не дороже 100 % цены по требованию
    jd = req["job_definition"]
    assert jd["timeout"]["attemptDurationSeconds"] == 12600
    assert jd["retryStrategy"]["evaluateOnExit"][0] == {"onStatusReason": "Host EC2*", "action": "RETRY"}
    assert jd["retryStrategy"]["evaluateOnExit"][-1] == {"onReason": "*", "action": "EXIT"}
    assert jd["containerProperties"]["logConfiguration"]["logDriver"] == "awslogs"
    assert jd["containerProperties"]["command"] == ["--job", "Ref::job_uri"]
    job = json.loads(files["s3/inputs/job.json"])
    assert job["time_limit_s"] == 12600 - 600 and job["np"] == 64
    assert job["run_prefix"] == "s3://PLACEHOLDER-bucket/PLACEHOLDER-prefix/gtu/gtu-EXAMPLE-kv3-char"
    assert req["submit_array_job"]["parameters"]["job_uri"] == job["run_prefix"] + "/inputs/job.json"


def test_manifest_records_code_solver_plan():
    import hashlib
    _, files, _ = render.render(cfg(), "us-east-1", "line", rates(), GIT)
    m = json.loads(files["s3/manifest.json"])
    assert m["code"]["git_sha"] == GIT["sha"]
    assert m["solver"]["version"] == "v8.6.0" and m["solver"]["commit"].startswith("7048fd35")
    assert m["plan"]["sha256"] == hashlib.sha256(files["s3/inputs/plan.csv"].encode()).hexdigest()
    assert m["plan"]["rows"] == 80 and m["plan"]["array_size"] == 10
    assert set(m["layout"]) == {"inputs/", "outputs/<row_id>/attempt-<n>/", "logs/<row_id>/attempt-<n>.log",
                                "manifest.json"}


def test_plan_index_to_row():
    rows = render.build_plan(cfg(), "line")
    assert len(rows) == 2 * 5 * 8
    assert [int(r["index"]) for r in rows] == list(range(80))
    assert len({r["row_id"] for r in rows}) == 80
    line3 = [r for r in rows if r["line"] == 3]
    assert len(line3) == 8 and {r["stage"] for r in line3} == {"kv3-c1"} and {r["speed_frac"] for r in line3} == {"1.00"}
    assert line3[0]["iter_scale"] == "1.0" and line3[1]["iter_scale"] == str(round(1000 / 6800, 4))
    assert float(line3[0]["p_out_frac"]) == pytest.approx(0.85) and float(line3[-1]["p_out_frac"]) == pytest.approx(1.12)
    assert all(r["example"] == "true" for r in rows)
    assert {r["iter_scale"] for r in render.build_plan(cfg(), "point")} == {"1.0"}


def _broken(mutate):
    c = cfg()
    req, _, _ = render.render(c, "us-east-1", "line", rates(), GIT)
    req = copy.deepcopy(req)
    mutate(req)
    return render.validate_api(req) + render.validate_limits(req, c, rates(), "us-east-1")


def test_validator_catches_bad_enum():
    errs = _broken(lambda r: r["compute_environment"]["computeResources"].__setitem__("allocationStrategy", "LOWEST"))
    assert any("allocationStrategy" in e for e in errs)


def test_validator_catches_missing_required_and_unknown_field():
    errs = _broken(lambda r: r["compute_environment"]["computeResources"].pop("maxvCpus"))
    assert any("maxvCpus" in e for e in errs)
    errs = _broken(lambda r: r["job_definition"].__setitem__("attemptDuration", 60))
    assert any("attemptDuration" in e for e in errs)


def test_validator_catches_documented_limits():
    errs = _broken(lambda r: r["compute_environment"]["computeResources"].__setitem__("maxvCpus", 384))
    assert any("квоты" in e for e in errs)
    errs = _broken(lambda r: r["job_definition"]["timeout"].__setitem__("attemptDurationSeconds", 30))
    assert any("60" in e for e in errs)
    errs = _broken(lambda r: r["job_definition"]["retryStrategy"].__setitem__("attempts", 11))
    assert any("attempts" in e for e in errs)
    errs = _broken(lambda r: r["submit_array_job"]["arrayProperties"].__setitem__("size", 1))
    assert any("size" in e for e in errs)


def test_render_main_rejects_other_region(tmp_path):
    assert render.main(["--region", "eu-west-2", "--out", str(tmp_path)]) == 1
    assert not any(tmp_path.iterdir())


# ================================================================ entrypoint: задание по JSON
def make_run(tmp_path, command, mode="line", time_limit=600, rows=None):
    root = tmp_path / "s3"
    prefix = "s3://bkt/pfx/gtu/run-test"
    base = root / "bkt/pfx/gtu/run-test"
    (base / "inputs").mkdir(parents=True)
    rows = rows or [{"index": i, "row_id": f"r{i}", "line": i // 3, "speed_frac": "0.90", "state": f"0.90,{i}"}
                    for i in range(6)]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    (base / "inputs/plan.csv").write_text(buf.getvalue())
    (base / "inputs/mesh.txt").write_text("mesh")
    job = {"schema": entrypoint.SCHEMA, "run_id": "run-test", "run_prefix": prefix, "plan": "inputs/plan.csv",
           "mode": mode, "command": command, "np": 2, "time_limit_s": time_limit}
    (base / "inputs/job.json").write_text(json.dumps(job))
    return root, prefix + "/inputs/job.json", base


def run_ep(root, job_uri, tmp_path, index, home=None):
    env_keep = os.environ.get("GTU_HOME")
    os.environ["GTU_HOME"] = str(home or tmp_path)
    try:
        return entrypoint.main(["--job", job_uri, "--index", str(index), "--local-root", str(root),
                                "--work", str(tmp_path / "work")])
    finally:
        if env_keep is None:
            os.environ.pop("GTU_HOME", None)
        else:
            os.environ["GTU_HOME"] = env_keep


CMD_OK = ('echo "row {row_id} f={speed_frac} np=$GTU_NP prev=$GTU_PREV_RESTART state=$GTU_PREV_STATE" '
          '> "$GTU_OUT/result.txt"; cp "$GTU_INPUTS/mesh.txt" "$GTU_OUT/"; touch "$GTU_OUT/restart_{row_id}.dat"')


def test_entrypoint_line_mode_uploads_each_point(tmp_path):
    root, uri, base = make_run(tmp_path, CMD_OK)
    assert run_ep(root, uri, tmp_path, 1) == 0
    for rid in ("r3", "r4", "r5"):
        st = json.loads((base / f"outputs/{rid}/attempt-1/_status.json").read_text())
        assert st["status"] == "ok" and st["exit_code"] == 0 and st["line"] == 1
        assert (base / f"outputs/{rid}/attempt-1/mesh.txt").read_text() == "mesh"
        assert (base / f"logs/{rid}/attempt-1.log").exists()
    r4 = (base / "outputs/r4/attempt-1/result.txt").read_text()
    assert "f=0.90 np=2" in r4 and "restart_r3.dat" in r4 and "state=0.90,3" in r4
    assert "prev= " in (base / "outputs/r3/attempt-1/result.txt").read_text()   # первая точка — без рестарта
    assert not (base / "outputs/r0").exists()


def test_entrypoint_line_stops_after_failure(tmp_path):
    cmd = 'if [ "{row_id}" = "r1" ]; then exit 7; fi; echo ok > "$GTU_OUT/x"'
    root, uri, base = make_run(tmp_path, cmd)
    assert run_ep(root, uri, tmp_path, 0) == 0          # линия дала точку r0 — partial, выход 0
    assert json.loads((base / "outputs/r1/attempt-1/_status.json").read_text())["exit_code"] == 7
    assert not (base / "outputs/r2").exists()            # после неудачи — skipped
    root2, uri2, _ = make_run(tmp_path / "b", 'exit 5')
    assert run_ep(root2, uri2, tmp_path / "b", 0) == 5   # первая точка не вышла — код команды


def test_entrypoint_point_mode_and_errors(tmp_path):
    root, uri, base = make_run(tmp_path, CMD_OK, mode="point")
    assert run_ep(root, uri, tmp_path, 4) == 0
    assert (base / "outputs/r4/attempt-1/_status.json").exists() and not (base / "outputs/r3").exists()
    assert run_ep(root, uri, tmp_path, 99) == entrypoint.EXIT_JOB_ERROR
    root2, uri2, _ = make_run(tmp_path / "c", "echo {no_such_field}")
    assert run_ep(root2, uri2, tmp_path / "c", 0) == entrypoint.EXIT_JOB_ERROR


def test_entrypoint_dry_run_writes_nothing(tmp_path, capsys):
    root, uri, base = make_run(tmp_path, CMD_OK)
    assert entrypoint.main(["--job", uri, "--index", "0", "--local-root", str(root), "--dry-run",
                            "--work", str(tmp_path / "w")]) == 0
    assert "r0: echo" in capsys.readouterr().out
    assert not (base / "outputs").exists()


def test_entrypoint_timeout_kills_silent_command(tmp_path):
    root, uri, base = make_run(tmp_path, "sleep 600", mode="point", time_limit=62)
    t0 = time.monotonic()
    rc = run_ep(root, uri, tmp_path, 0)
    assert time.monotonic() - t0 < 60
    st = json.loads((base / "outputs/r0/attempt-1/_status.json").read_text())
    assert st["status"] == "timeout" and rc != 0


def test_entrypoint_sigterm_uploads_and_exits_143(tmp_path):
    root, uri, base = make_run(tmp_path, 'echo started > "$GTU_OUT/s.txt"; sleep 600', mode="point")
    env = dict(os.environ, GTU_HOME=str(tmp_path))
    p = subprocess.Popen([sys.executable, str(HERE / "entrypoint.py"), "--job", uri, "--index", "0",
                          "--local-root", str(root), "--work", str(tmp_path / "work")], env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = tmp_path / "work/r0/out/s.txt"
    for _ in range(100):
        if out.exists():
            break
        time.sleep(0.1)
    p.send_signal(signal.SIGTERM)
    assert p.wait(timeout=30) == entrypoint.EXIT_TERMINATED
    st = json.loads((base / "outputs/r0/attempt-1/_status.json").read_text())
    assert st["status"] == "terminated" and (base / "outputs/r0/attempt-1/s.txt").exists()


def test_render_command_keeps_bash_expansions():
    t = 'run --f {speed_frac} ${GTU_PREV_RESTART:+--restart-from "$GTU_PREV_RESTART"} "$GTU_OUT"'
    assert entrypoint.render_command(t, {"speed_frac": "0.7"}) == \
        'run --f 0.7 ${GTU_PREV_RESTART:+--restart-from "$GTU_PREV_RESTART"} "$GTU_OUT"'
    with pytest.raises(KeyError):
        entrypoint.render_command("{x}", {})


def test_rendered_job_runs_through_entrypoint_dry(tmp_path):
    """job.json и plan.csv из render.py читаются точкой входа: 10 линий × 8 точек, команда собирается."""
    _, files, _ = render.render(cfg(), "us-east-1", "line", rates(), GIT)
    base = tmp_path / "s3/PLACEHOLDER-bucket/PLACEHOLDER-prefix/gtu/gtu-EXAMPLE-kv3-char"
    for rel in ("inputs/plan.csv", "inputs/job.json"):
        (base / rel).parent.mkdir(parents=True, exist_ok=True)
        (base / rel).write_text(files["s3/" + rel])
    job = entrypoint.load_job(json.loads(files["s3/inputs/job.json"])["run_prefix"] + "/inputs/job.json",
                              entrypoint.LocalStorage(tmp_path / "s3"))
    rows = entrypoint.select_rows(entrypoint.read_plan(files["s3/inputs/plan.csv"]), "line", 9)
    assert [r["row_id"] for r in rows][0] == "kv3-c2-n1.05-p1" and len(rows) == 8
    cmd = entrypoint.render_command(job["command"], rows[1])
    assert f"--iter-scale {round(1000 / 6800, 4)}" in cmd and "--frac 1.05" in cmd and "--p-out" in cmd


# ================================================================ смета
def load_estimate():
    c = cost_estimate.load_toml(HERE / "cost_estimate.toml")
    return c, cost_estimate.estimate(c, cost_estimate.load_toml(HERE / "batch.toml"),
                                     cost_estimate.load_toml(HERE / "aws_rates.toml"))


def test_results_json_is_current():
    """cost_estimate_results.json и CSV в репозитории совпадают с пересчётом (кроме строки версии Python)."""
    _, res = load_estimate()
    stored = json.loads((HERE / "cost_estimate_results.json").read_text())
    res = json.loads(json.dumps(res))
    stored.pop("software")
    res.pop("software")
    assert stored == res
    with open(HERE / "cost_estimate_lines.csv", newline="") as f:
        n_csv = sum(1 for _ in csv.DictReader(f))
    assert n_csv == len(res["lines"]) + sum(len(v["lines"]) for v in res["variants"].values())


def test_estimate_arithmetic_and_ceiling():
    c, res = load_estimate()
    b = tomllib.load(open(HERE / "batch.toml", "rb"))
    r = tomllib.load(open(HERE / "aws_rates.toml", "rb"))["regions"]["us-east-1"]["ec2"]
    s = res["scenarios"]["series_char"]
    od = max(r[t]["on_demand_usd_h"] for t in b["compute"]["instance_types"])
    assert od == r["c7a.16xlarge"]["on_demand_usd_h"] == 3.28448
    bound_h = b["job"]["attempt_duration_s"] / 3600 + (5 + 10) / 60
    comp = [x for x in res["lines"] if x["scenario"] == "series_char" and x["item"] == "compute"]
    ub = next(x for x in comp if x["level"] == "upper_bound")
    assert ub["cost_usd"] == pytest.approx(10 * b["job"]["retry_attempts"] * bound_h * od, rel=1e-6)
    odl = next(x for x in comp if x["level"] == "on_demand")
    assert odl["cost_usd"] == pytest.approx(10 * s["expected_h_per_job"] * od, rel=1e-6)
    for lv in cost_estimate.LEVELS:
        assert res["plan_total"][lv] == pytest.approx(sum(x["cost_usd"] for x in res["lines"] if x["level"] == lv))
    pt = res["plan_total"]
    assert pt["spot_estimate"] < pt["on_demand"] < pt["upper_bound"]
    assert res["ceiling_check"]["upper_bound_ok"] is (pt["upper_bound"] <= 300.0) is True
    assert s["iterations_per_job"] == 6800 + 7 * 1000 and s["jobs"] == 10 and s["points"] == 80


def test_estimate_has_all_line_kinds_with_sources():
    _, res = load_estimate()
    items = {x["item"] for x in res["lines"]}
    for need in ("compute", "ebs_root", "public_ipv4", "s3_storage", "s3_put", "s3_get", "transfer_out_internet",
                 "cw_logs_ingest_us-east-1", "cw_logs_storage_us-east-1", "ecr_storage_us-east-1"):
        assert need in items
    assert all(x["rate_source"].startswith("https://") for x in res["lines"])
    eu = res["variants"]["series_char_eu-north-1"]["lines"]
    assert any(x["item"] == "transfer_eu-north-1_to_us-east-1" for x in eu)
    for k in ("method", "software", "date", "inputs"):
        assert k in res
    assert res["date"] == "2026-11-14" and res["inputs"]["example"] is True


def test_iterations_match_su2_design_plan():
    """Итерации плана design в batch.toml и в смете — те же, что PLANS["design"] в run_case.py узла su2."""
    import importlib.util  # noqa: PLC0415
    spec = importlib.util.spec_from_file_location("run_case", HERE.parent / "su2" / "run_case.py")
    rc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rc)
    n = sum(int(s["iter"]) for s in rc.PLANS["design"])
    assert cfg()["plan"]["iters_first"] == n
    c, _ = load_estimate()
    assert {s["iterations"] for s in c["scenario"] if s["path"] == "runner"} == {n}


def test_estimate_bench_matches_runner_jobs_file():
    """Срок и диск замера в смете — те же, что в файле задач раннера; файл задач принимает сам runner.py."""
    sys.path.insert(0, str(HERE.parent / "aws"))
    import runner  # noqa: PLC0415
    jobs = {j.name: j for j in runner.load_jobs(HERE / "bench_runner_jobs.json")}
    c, _ = load_estimate()
    sc = {s["name"]: s for s in c["scenario"]}
    for jname, sname in (("bench-c7a16", "bench_c7a16"), ("bench-c8i8", "bench_c8i8")):
        j, s = jobs[jname], sc[sname]
        assert j.instance_type == s["instance_type"] and j.types == [s["instance_type"]]
        assert j.max_min == round(s["cap_h"] * 60) and j.disk_gb == s["root_gb"]
    assert jobs["bench-c8i8"].vcpu == 32


# ================================================================ тарифы
def test_rates_file_has_sources_and_consistent_spot():
    r = tomllib.load(open(HERE / "aws_rates.toml", "rb"))
    for key in ("AmazonEC2_us-east-1", "AmazonS3_us-east-1", "AmazonCloudWatch_us-east-1", "AmazonVPC_us-east-1",
                "AWSDataTransfer_us-east-1", "AmazonECR_us-east-1", "AmazonEC2_eu-north-1", "spot_advisor"):
        assert r["sources"][key]["url"].startswith("https://") and r["sources"][key]["publication"]
    for region, reg in r["regions"].items():
        for t, e in reg["ec2"].items():
            if e["spot_savings_pct"] >= 0:
                assert e["spot_estimate_usd_h"] == pytest.approx(e["on_demand_usd_h"] * (1 - e["spot_savings_pct"] / 100),
                                                                 abs=1e-6)
    use1 = r["regions"]["us-east-1"]
    assert use1["ec2"]["c8i.8xlarge"]["vcpu"] == 32 and use1["ec2"]["c7a.16xlarge"]["memory_gib"] == 128.0


def test_fetch_rates_extracts_from_fixture():
    """Разбор строк Price List на выдуманном образце (FIXTURE): нужная строка выбирается, лишние отсеиваются."""
    base = {"TermType": "OnDemand", "Tenancy": "Shared", "Operating System": "Linux", "Pre Installed S/W": "NA",
            "CapacityStatus": "Used", "License Model": "No License required", "Unit": "Hrs", "vCPU": "64",
            "Memory": "128 GiB", "Physical Processor": "FIXTURE CPU", "SKU": "FIXTURE-SKU", "EffectiveDate": "2026-11-01",
            "usageType": "BoxUsage"}
    rows = []
    for t in fetch_rates.INSTANCE_TYPES:
        rows.append({**base, "Instance Type": t, "PricePerUnit": "1.5"})
        rows.append({**base, "Instance Type": t, "PricePerUnit": "9.9", "Tenancy": "Dedicated"})
    rows.append({**base, "Instance Type": "", "usageType": "EBS:VolumeUsage.gp3", "PricePerUnit": "0.08",
                 "Unit": "GB-Mo", "PriceDescription": "FIXTURE gp3"})
    inst, gp3 = fetch_rates.extract_ec2(rows)
    assert inst["c7a.16xlarge"]["on_demand_usd_h"] == 1.5 and gp3["usd"] == 0.08
    offer = {"products": {"A": {"attributes": {"usagetype": "X-PublicIPv4:InUseAddress"}},
                          "B": {"attributes": {"usagetype": "other"}}},
             "terms": {"OnDemand": {"A": {"o": {"priceDimensions": {"d": {"pricePerUnit": {"USD": "0.005"}, "unit": "Hrs",
                                                                         "beginRange": "0", "description": "FIXTURE"}}}}}}}
    assert fetch_rates.pick(offer, lambda a: a["usagetype"].endswith("PublicIPv4:InUseAddress"))["usd"] == 0.005


def test_dump_toml_roundtrip():
    data = {"meta": {"a": 1, "b": "строка"}, "regions": {"us-east-1": {"ec2": {"c7a.16xlarge": {"x": 1.5, "l": [1, 2]}}}}}
    assert tomllib.loads(fetch_rates.dump_toml(data)) == data


# ================================================================ Dockerfile (статически; сборка — README)
def test_dockerfile_static():
    d = (HERE.parent / "docker/Dockerfile").read_text()
    b = tomllib.load(open(HERE / "batch.toml", "rb"))
    assert d.count("\nFROM ") == 2 and "AS su2-build" in d and "AS runtime" in d
    assert "ARG SU2_TAG=v8.6.0" in d and f"ARG SU2_COMMIT={b['solver']['commit']}" in d
    assert '"/opt/gtu/routes/gtu/platform/batch/entrypoint.py"' in d and (HERE / "entrypoint.py").exists()
    assert "calculix-ccx" in d and "-Dwith-mpi=enabled" in d and "USER 1000:1000" in d
    assert "OMPI_MCA_osc=^ucx" in d   # без него SU2_CFD под mpirun в контейнере падает в MPI_Win_create
    req = (HERE.parent / "docker/requirements.txt").read_text().split("\n")
    pkgs = [x for x in req if x and not x.startswith("#")]
    assert {p.split("==")[0] for p in pkgs} == {"cantera", "numpy", "scipy", "gmsh", "boto3"}
    assert all("==" in p for p in pkgs)
    ign = (HERE.parent / "docker/Dockerfile.dockerignore").read_text()
    assert "!routes/gtu/" in ign and "routes/gtu/**/runs/" in ign
