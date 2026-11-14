"""Тесты runner.py без сети: botocore Stubber вместо AWS, user-data — на подменённых командах.

  python3 -m pytest -q routes/gtu/platform/aws/test_runner.py

Ключи в тестах поддельные (ставятся в окружение на время теста); любой вызов AWS,
для которого нет заготовленного ответа, Stubber превращает в ошибку — сети нет.
"""
import datetime as dt
import io
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import botocore.session
import pytest
from botocore.response import StreamingBody
from botocore.stub import ANY, Stubber

sys.path.insert(0, str(Path(__file__).resolve().parent))
import runner as R  # noqa: E402

FAKE_AKID = "AKIASYNTHETIC0000000"  # выдуманный ключ: пометка SYNTHETIC по правилу gate.yml
FAKE_SECRET = "fakeSecretValueForTestsOnly/0123456789abcdEF"
BUCKET = "alpha-o-results-123456789012"
PREFIX = "alpha-o"
RUN = "20261114T000000Z-ab12"
T0 = 1_800_000_000.0
UTC = dt.timezone.utc


@pytest.fixture(autouse=True)
def fake_env(monkeypatch):
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", FAKE_AKID)
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", FAKE_SECRET)
    monkeypatch.delenv("AWS_SESSION_TOKEN", raising=False)
    monkeypatch.setenv("AWS_CONFIG_FILE", "/nonexistent")
    monkeypatch.setenv("AWS_SHARED_CREDENTIALS_FILE", "/nonexistent")
    monkeypatch.setenv("AWS_EC2_METADATA_DISABLED", "true")


class Clock:
    def __init__(self, t=T0):
        self.t = t

    def __call__(self):
        return self.t

    def sleep(self, s):
        self.t += s


def client(name):
    return botocore.session.get_session().create_client(
        name, region_name="us-east-1", aws_access_key_id="testing", aws_secret_access_key="testing")


class Aws:
    """Клиенты со Stubber; в конце теста — проверка, что все заготовки израсходованы."""

    def __init__(self):
        self.ec2, self.s3, self.ssm = client("ec2"), client("s3"), client("ssm")
        self.st = {n: Stubber(getattr(self, n)) for n in ("ec2", "s3", "ssm")}
        for s in self.st.values():
            s.activate()

    def done(self):
        for s in self.st.values():
            s.assert_no_pending_responses()


def make_cfg(tmp_path, **kw):
    base = dict(bucket=BUCKET, prefix=PREFIX, instance_profile="gtu-test-profile", run_id=RUN,
                ledger_path=tmp_path / "ledger.jsonl", out_dir=tmp_path / "runs", poll_s=60.0)
    base.update(kw)
    return R.Config(**base)


def make_job(**kw):
    base = dict(name="su2-comp130", command='echo hi > "$GTU_OUT/a.txt"', max_min=60)
    base.update(kw)
    job = R.Job(**base)
    R.validate_job(job)
    return job


def spot_params(itype, now=T0, window_h=6.0):
    end = dt.datetime.fromtimestamp(now, tz=UTC)
    return {"InstanceTypes": [itype], "ProductDescriptions": ["Linux/UNIX"],
            "StartTime": end - dt.timedelta(hours=window_h), "EndTime": end}


def spot_resp(itype, prices):
    items = [{"AvailabilityZone": az, "InstanceType": itype, "ProductDescription": "Linux/UNIX",
              "SpotPrice": f"{p:.4f}", "Timestamp": dt.datetime.fromtimestamp(T0 + off, tz=UTC)}
             for az, lst in prices.items() for p, off in lst]
    return {"SpotPriceHistory": items}


def body(b):
    return StreamingBody(io.BytesIO(b), len(b))


def write_ledger(path, recs):
    path.write_text("".join(json.dumps(r) + "\n" for r in recs), encoding="utf-8")


def inst(iid, state="running", tags=None, az="us-east-1b", reason=None):
    d = {"InstanceId": iid, "State": {"Name": state}, "Placement": {"AvailabilityZone": az}}
    if tags is not None:
        d["Tags"] = [{"Key": k, "Value": v} for k, v in tags.items()]
    if reason:
        d["StateReason"] = {"Code": reason, "Message": reason}
    return d


def own_tags(job="su2-comp130", run=RUN):
    return {"Project": R.DEFAULT_TAG, "RunId": run, "Job": job, "Name": f"gtu-{job}"}


def reap_filters(run=RUN):
    return {"Filters": [{"Name": "tag:Project", "Values": [R.DEFAULT_TAG]},
                        {"Name": "tag:RunId", "Values": [run]},
                        {"Name": "instance-state-name", "Values": R.LIVE_STATES}]}


PRICES_C7A = {"us-east-1a": [(0.9, -3000), (1.1, -100)], "us-east-1b": [(1.3, -2000), (1.2, -50)]}
HIST = [
    # закрытая: (1.00 + 0.01) долл./ч × 1 ч = 1.01
    {"ev": "launch", "run_id": "old", "job": "x", "iid": "i-closed", "itype": "c7a.16xlarge", "az": "us-east-1a",
     "spot_usd_h": 1.0, "max_usd_h": 1.5, "ebs_usd_h": 0.01, "t_launch": T0 - 7200, "allot_min": 70},
    {"ev": "end", "run_id": "old", "job": "x", "iid": "i-closed", "t_end": T0 - 3600},
    # открытая: верхняя граница (2.00 + 0.01) × 70/60 = 2.345
    {"ev": "launch", "run_id": "old", "job": "y", "iid": "i-open", "itype": "c7a.16xlarge", "az": "us-east-1a",
     "spot_usd_h": 1.0, "max_usd_h": 2.0, "ebs_usd_h": 0.01, "t_launch": T0 - 600, "allot_min": 70},
]


# ----------------------------------------------------------------------------- бюджет
def test_budget_estimate_numbers(tmp_path):
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               spot_params("c7a.16xlarge"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    job = make_job()
    info = r.price("c7a.16xlarge")
    assert info.latest == {"us-east-1a": 1.1, "us-east-1b": 1.2} and info.window_max == 1.3
    ok, maxp, nums = r.budget_check(job, "c7a.16xlarge", info)
    assert maxp == pytest.approx(1.56)                      # ceil(1.3 × 1.2, 4 знака)
    ebs_h = 50 * 0.08 / 730
    new = (1.56 + ebs_h + 0.005) * 70 / 60                  # + публичный IPv4; max_min 60 + запас 10
    assert ok
    assert nums["spent"] == pytest.approx(1.01, abs=1e-4)
    assert nums["committed"] == pytest.approx(2.345, abs=1e-4)
    assert nums["new"] == pytest.approx(new, abs=1e-4)
    assert nums["total"] == pytest.approx(1.01 + 2.345 + new, abs=1e-3)
    rec = [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "budget"][-1]
    assert rec["ok"] is True and rec["max_usd_h"] == pytest.approx(1.56) and rec["job"] == job.name
    aws.done()


def test_budget_refusal_blocks_launch(tmp_path):
    # Израсходовано 299 долл. (закрытая машина 2.99 долл./ч × 100 ч); новая задача ≈ 1.8 → отказ.
    write_ledger(tmp_path / "ledger.jsonl", [
        {"ev": "launch", "run_id": "old", "job": "x", "iid": "i-big", "itype": "c7a.16xlarge", "az": "us-east-1a",
         "spot_usd_h": 2.99, "max_usd_h": 3.0, "ebs_usd_h": 0.0, "t_launch": T0 - 400 * 3600, "allot_min": 6000},
        {"ev": "end", "run_id": "old", "job": "x", "iid": "i-big", "t_end": T0 - 300 * 3600}])
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", {"us-east-1a": [(1.4, -60)]}),
                   spot_params("c7i.16xlarge"))
    e.add_response("describe_instances", {"Reservations": []}, reap_filters())   # сборщик в finally
    # run_instances НЕ заготовлен: если исполнитель попробует запустить машину, Stubber упадёт.
    # потолок задан явно: тест проверяет отказ при превышении, а не значение по умолчанию
    r = R.Runner(make_cfg(tmp_path, ceiling=300.0), aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res == {"su2-comp130": {"status": "refused", "attempts": 1}}
    budget = [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "budget"]
    assert [(b["itype"], b["ok"]) for b in budget] == [("c7a.16xlarge", False), ("c7i.16xlarge", False)]
    assert all(b["total"] > b["ceiling"] for b in budget)
    assert not [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "launch" and x["run_id"] == RUN]
    aws.done()


def test_ceiling_counts_open_instances_and_parallel_launches(tmp_path):
    """Открытые машины входят в прогноз целиком — параллельный запуск не обходит потолок."""
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    s = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s["open"] == ["i-open"] and s["closed"] == 1
    assert s["spent"] == pytest.approx(1.01) and s["committed"] == pytest.approx(2.345)


def test_second_parallel_launch_refused_by_first_commitment(tmp_path):
    """Потолок 3 долл.: первая задача (≈1.83) запущена, вторая с её верхней границей — отказ."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_response("run_instances", run_resp("i-a"), expected_run_params(job="a"))
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", {"us-east-1a": [(1.4, -60)]}),
                   spot_params("c7i.16xlarge"))
    r = R.Runner(make_cfg(tmp_path, ceiling=3.0), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    from collections import deque
    pending = deque([(make_job(name="a"), 1), (make_job(name="b"), 1)])
    r._launch_some(pending)
    assert list(r.active) == ["i-a"] and r.results["b"]["status"] == "refused"
    budget = [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "budget"]
    assert [(b["job"], b["ok"]) for b in budget] == [("a", True), ("b", False), ("b", False)]
    assert budget[1]["committed"] == pytest.approx(budget[0]["new"], abs=1e-3)
    aws.done()


# ----------------------------------------------------------------------------- user-data
def test_userdata_watchdog_first_and_no_keys():
    job = make_job(max_min=90)
    ud = R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="claude/gtu-x",
                          bucket=BUCKET, key_prefix=f"{PREFIX}/gtu/{RUN}/{job.name}")
    lines = ud.splitlines()
    assert lines[0] == "#!/bin/bash"
    assert lines[1] == "shutdown -P +90"                       # сторож — первая команда
    commands = [x for x in lines[1:] if x.strip() and not x.lstrip().startswith("#")]
    assert commands[0] == "shutdown -P +90"
    assert FAKE_AKID not in ud and FAKE_SECRET not in ud
    low = ud.lower()
    for w in ("aws_secret_access_key", "aws_access_key_id", "aws_session_token", "akia"):
        assert w not in low
    assert "poweroff" in ud and "trap finish EXIT" in ud
    assert f"s3_key={PREFIX}/gtu/{RUN}/{job.name}".lower() in low.replace("'", "")
    assert len(ud.encode()) <= R.USERDATA_LIMIT


@pytest.mark.parametrize("image", ["ubuntu2404", "al2023"])
def test_userdata_bash_syntax(tmp_path, image):
    kw = {"apt": ["calculix-ccx", "gmsh"]} if image == "ubuntu2404" else {"dnf": ["gcc"]}
    job = make_job(image=image, pip=["numpy>=1.26", "gmsh"], **kw)
    ud = R.build_userdata(job, run_id=RUN, attempt=2, repo=R.DEFAULT_REPO, branch="main",
                          bucket=BUCKET, key_prefix=f"{PREFIX}/gtu/{RUN}/{job.name}")
    f = tmp_path / "ud.sh"
    f.write_text(ud)
    assert subprocess.run(["bash", "-n", str(f)]).returncode == 0
    assert "'numpy>=1.26'" in ud                                # спецификации пакетов закавычены


def test_secrets_refused():
    with pytest.raises(R.RunnerError):
        make_job(command=f"export AWS_ACCESS_KEY_ID={FAKE_AKID}; run")
    with pytest.raises(R.RunnerError):
        make_job(command=f"echo {FAKE_SECRET}")
    with pytest.raises(R.RunnerError):
        make_job(command="echo AKIAEXAMPLE000000000")
    job = R.Job(name="j", command="echo ok", max_min=60)
    R.validate_job(job)
    job.command = f"echo {FAKE_SECRET}"                         # в обход проверки задачи
    with pytest.raises(R.RunnerError):
        R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main",
                         bucket=BUCKET, key_prefix="p/gtu/r/j")
    with pytest.raises(R.RunnerError):                          # предел 16 КБ
        R.build_userdata(make_job(command="echo " + "x" * 20000), run_id=RUN, attempt=1,
                         repo=R.DEFAULT_REPO, branch="main", bucket=BUCKET, key_prefix="p/gtu/r/j")
    with pytest.raises(R.RunnerError):
        R.build_userdata(make_job(), run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main; rm -rf /",
                         bucket=BUCKET, key_prefix="p/gtu/r/j")


SHIMS = {
    "shutdown": 'echo "shutdown $*" >> "$CALLS"\n',
    "poweroff": 'echo poweroff >> "$CALLS"\n',
    "timeout": 'echo "timeout $*" >> "$CALLS"\nexec "$REALTIMEOUT" "$@"\n',
    "apt-get": 'echo "apt-get $*" >> "$CALLS"\ncase " $* " in *" failpkg "*) exit 100;; esac\n',
    "git": 'echo "git $*" >> "$CALLS"\nfor last; do :; done\n[ "$1" = clone ] && mkdir -p "$last"\nexit 0\n',
    "python3": ('if [ "$1" = -m ] && [ "$2" = venv ]; then\n'
                '  mkdir -p "$3/bin"\n'
                "  printf '#!/bin/bash\\nexec env PYTHONPATH=\"$FAKEMODS\" \"$REALPY\" \"$@\"\\n' > \"$3/bin/python\"\n"
                "  printf '#!/bin/bash\\necho \"pip $*\" >> \"$CALLS\"\\n' > \"$3/bin/pip\"\n"
                '  chmod +x "$3/bin/python" "$3/bin/pip"; echo venv >> "$CALLS"; exit 0\n'
                'fi\nexec "$REALPY" "$@"\n'),
}
FAKE_BOTO3 = '''import os, shutil
class _S3:
    def upload_file(self, path, bucket, key):
        dst = os.path.join(os.environ["S3DIR"], bucket, key)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(path, dst)
        with open(os.environ["CALLS"], "a") as fh:
            fh.write("upload " + key + "\\n")
def client(name, *a, **k):
    assert name == "s3"
    return _S3()
'''


def run_box(tmp_path, job):
    """Выполнить user-data локально: shutdown/poweroff/apt-get/git/python3 подменены, boto3 — поддельный."""
    shims, mods, s3dir, calls = tmp_path / "shims", tmp_path / "mods" / "boto3", tmp_path / "s3", tmp_path / "calls"
    shims.mkdir()
    mods.mkdir(parents=True)
    for name, text in SHIMS.items():
        (shims / name).write_text("#!/bin/bash\n" + text)
        (shims / name).chmod(0o755)
    (mods / "__init__.py").write_text(FAKE_BOTO3)
    kp = f"pfx/gtu/rid1/{job.name}"
    ud = R.build_userdata(job, run_id="rid1", attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket="bkt",
                          key_prefix=kp, _root=str(tmp_path / "box"), _log=str(tmp_path / "box.log"))
    (tmp_path / "ud.sh").write_text(ud)
    env = {"PATH": f"{shims}:{os.environ['PATH']}", "CALLS": str(calls), "S3DIR": str(s3dir),
           "FAKEMODS": str(tmp_path / "mods"), "REALPY": sys.executable, "HOME": str(tmp_path),
           "REALTIMEOUT": shutil.which("timeout")}
    subprocess.run(["bash", str(tmp_path / "ud.sh")], env=env, timeout=120, check=False)
    return calls.read_text().splitlines(), s3dir / "bkt" / kp


def test_userdata_runs_locally_success(tmp_path):
    job = make_job(name="box-ok", max_min=30, pip=["cantera"],
                   command='echo hello > "$GTU_OUT/result.txt"\nmkdir -p "$GTU_OUT/sub"\n'
                           'echo 42 > "$GTU_OUT/sub/x.txt"\npwd > "$GTU_OUT/pwd.txt"')
    calls, s3 = run_box(tmp_path, job)
    assert calls[0] == "shutdown -P +30"                       # сторож первым
    assert calls[-1] == "poweroff"                              # гасится сама
    assert calls[-2] == "upload pfx/gtu/rid1/box-ok/_DONE.json" # метка — последней
    assert any(c.startswith("git clone -q --depth 1 -b main https://github.com/") for c in calls)
    assert "pip install -q cantera" in calls
    assert (s3 / "result.txt").read_text() == "hello\n"
    assert (s3 / "sub" / "x.txt").read_text() == "42\n"
    assert (s3 / "pwd.txt").read_text().strip() == str(tmp_path / "box" / "repo")
    assert (s3 / "_job.log").exists() and (s3 / "_stdout.log").exists()
    done = json.loads((s3 / "_DONE.json").read_text())
    assert done["phase"] == "done" and done["exit_code"] == 0 and done["run_id"] == "rid1"
    # выгрузка укладывается до сторожа: результаты — не дольше (срок − прошло − 120 с), метка — 60 с
    ups = [c.split()[1] for c in calls if c.startswith("timeout ") and "upload.py" in c]
    assert len(ups) == 2 and ups[1] == "60s"
    assert 30 * 60 - 120 - 60 <= int(ups[0].rstrip("s")) <= 30 * 60 - 120


def test_upload_window_fits_inside_watchdog():
    """Худший случай: команда снята по timeout -k 60 в последний момент — метка всё равно до сторожа."""
    for max_min in (20, 60, 720):
        job = make_job(max_min=max_min)
        ud = R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main",
                              bucket=BUCKET, key_prefix="p/gtu/r/j")
        reserve = 8 * 60
        cmd_end = max_min * 60 - reserve + 60                   # DEADLINE + -k 60 (от T0)
        left = min(1800, max(30, max_min * 60 - cmd_end - 120))  # как в finish()
        assert cmd_end + left + 60 <= max_min * 60 - 60         # результаты + метка, минута до сторожа
        assert 'up "$OUT" "$LEFT"' in ud and 'up "$G/_DONE.json" 60' in ud
        assert f"WD_END=$((T0 + {max_min} * 60))" in ud


def test_userdata_runs_locally_failures(tmp_path):
    calls, s3 = run_box(tmp_path / "a", _mk(tmp_path / "a", make_job(name="box-fail", max_min=30, command="exit 3")))
    done = json.loads((s3 / "_DONE.json").read_text())
    assert (done["phase"], done["exit_code"]) == ("done", 3)
    assert calls[0] == "shutdown -P +30" and calls[-1] == "poweroff"
    # сбой пакетов задачи: журнал и метка всё равно выгружены, клонирования нет
    calls, s3 = run_box(tmp_path / "b", _mk(tmp_path / "b", make_job(name="box-pkg", max_min=30, apt=["failpkg"])))
    done = json.loads((s3 / "_DONE.json").read_text())
    assert (done["phase"], done["exit_code"]) == ("packages", 16)
    assert not any(c.startswith("git clone") for c in calls)
    assert (s3 / "_job.log").exists() and calls[-1] == "poweroff"


def _mk(path, job):
    path.mkdir(parents=True)
    return job


# ----------------------------------------------------------------------------- run_instances
def expected_run_params(itype="c7a.16xlarge", maxp="1.5600", job="su2-comp130", disk=50, root="/dev/sda1"):
    tags = [{"Key": "Project", "Value": R.DEFAULT_TAG}, {"Key": "RunId", "Value": RUN},
            {"Key": "Job", "Value": job}, {"Key": "Name", "Value": f"gtu-{job}"}]
    return {
        "ImageId": "ami-0123456789abcdef0", "InstanceType": itype, "MinCount": 1, "MaxCount": 1,
        "IamInstanceProfile": {"Name": "gtu-test-profile"},
        "InstanceInitiatedShutdownBehavior": "terminate",
        "InstanceMarketOptions": {"MarketType": "spot",
                                  "SpotOptions": {"MaxPrice": maxp, "SpotInstanceType": "one-time",
                                                  "InstanceInterruptionBehavior": "terminate"}},
        "MetadataOptions": {"HttpTokens": "required", "HttpPutResponseHopLimit": 1, "HttpEndpoint": "enabled"},
        "BlockDeviceMappings": [{"DeviceName": root, "Ebs": {"VolumeSize": disk, "VolumeType": "gp3",
                                                             "DeleteOnTermination": True}}],
        "TagSpecifications": [{"ResourceType": "instance", "Tags": tags}, {"ResourceType": "volume", "Tags": tags}],
        "UserData": ANY, "ClientToken": ANY,
    }


def add_ami(aws, image="ubuntu2404"):
    aws.st["ssm"].add_response("get_parameter",
                               {"Parameter": {"Name": R.IMAGES[image]["ssm"], "Type": "String",
                                              "Value": "ami-0123456789abcdef0", "Version": 1}},
                               {"Name": R.IMAGES[image]["ssm"]})


def run_resp(iid, itype="c7a.16xlarge", az="us-east-1b"):
    return {"Instances": [{"InstanceId": iid, "InstanceType": itype, "Placement": {"AvailabilityZone": az},
                           "State": {"Name": "pending"}}]}


def test_run_instances_params(tmp_path):
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               spot_params("c7a.16xlarge"))
    add_ami(aws)
    aws.st["ec2"].add_response("run_instances", run_resp("i-0aaa"), expected_run_params())
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and act.iid == "i-0aaa" and act.az == "us-east-1b"
    assert act.deadline == T0 + 70 * 60
    p = r.last_params
    assert p["InstanceMarketOptions"]["MarketType"] == "spot"
    assert p["MetadataOptions"]["HttpTokens"] == "required"
    assert p["InstanceInitiatedShutdownBehavior"] == "terminate"
    assert p["UserData"].splitlines()[1] == "shutdown -P +60"
    assert len(p["ClientToken"]) <= 64 and "KeyName" not in p
    rec = [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "launch"][0]
    assert rec["iid"] == "i-0aaa" and rec["spot_usd_h"] == 1.2 and rec["max_usd_h"] == pytest.approx(1.56)
    assert rec["allot_min"] == 70
    aws.done()


def test_fallback_type_on_capacity_error(tmp_path):
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_client_error("run_instances", service_error_code="InsufficientInstanceCapacity",
                       service_message="no capacity", http_status_code=500, expected_params=expected_run_params())
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", {"us-east-1c": [(1.0, -60)]}),
                   spot_params("c7i.16xlarge"))
    e.add_response("run_instances", run_resp("i-0bbb", "c7i.16xlarge", "us-east-1c"),
                   expected_run_params(itype="c7i.16xlarge", maxp="1.2000"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and act.itype == "c7i.16xlarge"
    evs = [(x["ev"], x.get("itype"), x.get("code")) for x in R.Ledger(tmp_path / "ledger.jsonl").records()]
    assert ("launch_error", "c7a.16xlarge", "InsufficientInstanceCapacity") in evs
    assert ("launch", "c7i.16xlarge", None) in evs
    aws.done()


def test_dryrun_does_not_launch(tmp_path):
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               spot_params("c7a.16xlarge"))
    add_ami(aws)
    exp = dict(expected_run_params(), DryRun=True)
    aws.st["ec2"].add_client_error("run_instances", service_error_code="DryRunOperation",
                                   service_message="Request would have succeeded", http_status_code=412,
                                   expected_params=exp)
    cfg = make_cfg(tmp_path, dryrun=True)
    r = R.Runner(cfg, aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res["su2-comp130"]["status"] == "dryrun_ok"
    assert not (tmp_path / "ledger.jsonl").exists()             # общий журнал не тронут
    side = R.Ledger(cfg.out_dir / RUN / "journal.dryrun.jsonl").records()
    assert any(x["ev"] == "dryrun" and x["ok"] for x in side)
    aws.done()                                                  # ни terminate, ни сборщика не было


# ----------------------------------------------------------------------------- сборщик
def test_reaper_terminates_only_own(tmp_path):
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-own1", tags=own_tags()),
        inst("i-sibling", tags=own_tags(run="20261113T000000Z-ffff")),     # соседний прогон, тот же Project
        inst("i-untagged"),
        inst("i-foreign", tags={"Project": "other", "RunId": RUN}),
    ]}]}, reap_filters())
    aws.st["ec2"].add_response("terminate_instances", {"TerminatingInstances": []},
                               {"InstanceIds": ["i-own1", "i-own2"]})
    out = []
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=out.append)
    got = r.reap("тест", known=["i-own2"])
    assert got == ["i-own1", "i-own2"]
    assert any("i-sibling" in m and "не трогаю" in m for m in out)
    aws.done()


def test_reaper_retries_one_by_one_when_batch_rejected(tmp_path):
    """Один негодный id (свежая машина, NotFound) не оставляет живыми остальные; её запись — открыта."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-own1", tags=own_tags()), inst("i-own2", tags=own_tags())]}]}, reap_filters())
    e.add_client_error("terminate_instances", service_error_code="InvalidInstanceID.NotFound", http_status_code=400,
                       expected_params={"InstanceIds": ["i-fresh", "i-own1", "i-own2"]})
    e.add_client_error("terminate_instances", service_error_code="InvalidInstanceID.NotFound", http_status_code=400,
                       expected_params={"InstanceIds": ["i-fresh"]})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-own1"]})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-own2"]})
    out = []
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=out.append)
    job = make_job()
    for iid in ("i-fresh", "i-own1"):
        r.active[iid] = R.Active(job, 1, iid, "c7a.16xlarge", "us-east-1a", T0, T0 + 4200)
        r.ledger.append("launch", run_id=RUN, job=job.name, attempt=1, iid=iid, itype="c7a.16xlarge",
                        az="us-east-1a", spot_usd_h=1.0, max_usd_h=1.5, ebs_usd_h=0.0, ip_usd_h=0.0,
                        t_launch=T0, allot_min=70)
    r.reap("тест", known=list(r.active))                        # как в finally прогона
    assert any("погасите вручную" in m and "i-fresh" in m for m in out)
    s = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s["open"] == ["i-fresh"] and s["closed"] == 1        # неподтверждённое гашение не закрыто
    aws.done()


def test_cli_reap_closes_ledger_records(tmp_path, monkeypatch):
    monkeypatch.setenv("ALPHA_O_S3", f"s3://{BUCKET}/{PREFIX}")
    monkeypatch.setenv("ALPHA_O_INSTANCE_PROFILE", "gtu-test-profile")
    led = tmp_path / "ledger.jsonl"
    write_ledger(led, [dict(HIST[2], run_id=RUN)])
    aws = Aws()
    e = aws.st["ec2"]
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-open", tags=own_tags())]}]},
                   reap_filters())
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-open"]})
    e.add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-open", "shutting-down", own_tags())]}]}, {"Filters": [{"Name": "instance-id", "Values": ["i-open"]}]})
    monkeypatch.setattr(R, "_clients", lambda region: (aws.ec2, aws.s3, aws.ssm, None))
    assert R.main(["reap", "--run-id", RUN, "--ledger", str(led)]) == 0
    s = R.Ledger(led).summary()
    assert s["open"] == [] and s["closed"] == 1
    aws.done()


def test_client_token_differs_between_processes_with_same_run_id(tmp_path):
    job = make_job()
    toks = {R.Runner(make_cfg(tmp_path), out=lambda m: None).build_params(job, "c7a.16xlarge", "ami-x", 1.5, 1)
            ["ClientToken"] for _ in range(2)}
    assert len(toks) == 2 and all(len(t) <= 64 for t in toks)


def test_branch_check_exact_ref(monkeypatch):
    class P:
        returncode = 0
        stdout = "0123abcd\trefs/heads/xy/x\n"
    monkeypatch.setattr(R.subprocess, "run", lambda *a, **k: P())
    with pytest.raises(R.RunnerError):                          # ветки x нет, есть xy/x
        R._check_branch(R.DEFAULT_REPO, "x")
    R._check_branch(R.DEFAULT_REPO, "xy/x")


def test_collector_downloads_and_rejects_traversal(tmp_path):
    aws, clk = Aws(), Clock()
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    aws.st["s3"].add_response("list_objects_v2", {"IsTruncated": False, "Contents": [
        {"Key": pfx + "a.json"}, {"Key": pfx + "../evil.sh"}, {"Key": pfx + "sub/b.txt"}]},
        {"Bucket": BUCKET, "Prefix": pfx})
    aws.st["s3"].add_response("get_object", {"Body": body(b'{"a": 1}')}, {"Bucket": BUCKET, "Key": pfx + "a.json"})
    aws.st["s3"].add_response("get_object", {"Body": body(b"bbb")}, {"Bucket": BUCKET, "Key": pfx + "sub/b.txt"})
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    got = r.collect("su2-comp130")
    d = tmp_path / "runs" / RUN / "su2-comp130"
    assert got == ["a.json", "sub/b.txt"]
    assert (d / "a.json").read_text() == '{"a": 1}' and (d / "sub" / "b.txt").read_text() == "bbb"
    assert not (tmp_path / "runs" / RUN / "evil.sh").exists()
    aws.done()


# ----------------------------------------------------------------------------- полный цикл
def test_full_cycle_launch_poll_collect_reap_no_secrets(tmp_path):
    aws, clk, out = Aws(), Clock(), []
    e, s = aws.st["ec2"], aws.st["s3"]
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_response("run_instances", run_resp("i-0ccc"), expected_run_params())
    filt = {"Filters": [{"Name": "instance-id", "Values": ["i-0ccc"]}]}
    # опрос 1: машина работает, метки нет
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-0ccc", tags=own_tags())]}]}, filt)
    s.add_client_error("get_object", service_error_code="NoSuchKey", http_status_code=404,
                       expected_params={"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    # опрос 2: метка есть -> сбор -> гашение
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-0ccc", tags=own_tags())]}]}, filt)
    done = json.dumps({"run_id": RUN, "job": "su2-comp130", "phase": "done", "exit_code": 0}).encode()
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    s.add_response("list_objects_v2", {"IsTruncated": False, "Contents": [
        {"Key": pfx + "a.txt"}, {"Key": pfx + "_DONE.json"}]}, {"Bucket": BUCKET, "Prefix": pfx})
    s.add_response("get_object", {"Body": body(b"hi\n")}, {"Bucket": BUCKET, "Key": pfx + "a.txt"})
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-0ccc"]})
    # finally: сборщик по тегам — живых своих нет, чужая не трогается
    e.add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-sibling", tags=own_tags(run="other-run"))]}]}, reap_filters())
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=out.append)
    res = r.run([make_job()])
    assert res["su2-comp130"]["status"] == "ok" and res["su2-comp130"]["files"] == 2
    assert (tmp_path / "runs" / RUN / "su2-comp130" / "a.txt").read_text() == "hi\n"
    s_ = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s_["open"] == [] and s_["closed"] == 1
    assert s_["spent"] == pytest.approx((1.2 + 50 * 0.08 / 730 + 0.005) * 120 / 3600, abs=1e-4)   # 2 опроса по 60 с
    aws.done()
    # ключи и имя бакета (в нём номер аккаунта) нигде не напечатаны и не записаны
    texts = ["\n".join(out), (tmp_path / "ledger.jsonl").read_text(),
             (tmp_path / "runs" / RUN / "summary.json").read_text()]
    for t in texts:
        assert FAKE_AKID not in t and FAKE_SECRET not in t and BUCKET not in t and "123456789012" not in t


def test_deadline_reaps_hung_instance(tmp_path):
    aws, clk = Aws(), Clock()
    e, s = aws.st["ec2"], aws.st["s3"]
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_response("run_instances", run_resp("i-0ddd"), expected_run_params())
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-0ddd", tags=own_tags())]}]},
                   {"Filters": [{"Name": "instance-id", "Values": ["i-0ddd"]}]})
    s.add_client_error("get_object", service_error_code="NoSuchKey", http_status_code=404,
                       expected_params={"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-0ddd"]})
    s.add_response("list_objects_v2", {"IsTruncated": False}, {"Bucket": BUCKET, "Prefix": pfx})
    e.add_response("describe_instances", {"Reservations": []}, reap_filters())
    cfg = make_cfg(tmp_path, poll_s=71 * 60)                    # первый опрос — уже после срока 70 мин
    r = R.Runner(cfg, aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res["su2-comp130"]["status"] == "deadline"
    aws.done()


def test_lost_instance_is_retried_once(tmp_path):
    """Машина погашена без метки (прерывание спота) -> одна повторная попытка."""
    aws, clk = Aws(), Clock()
    e, s = aws.st["ec2"], aws.st["s3"]
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_response("run_instances", run_resp("i-1"), expected_run_params())
    e.add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-1", "terminated", own_tags(), reason="Server.SpotInstanceTermination")]}]},
        {"Filters": [{"Name": "instance-id", "Values": ["i-1"]}]})
    s.add_client_error("get_object", service_error_code="NoSuchKey", http_status_code=404,
                       expected_params={"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    s.add_response("list_objects_v2", {"IsTruncated": False}, {"Bucket": BUCKET, "Prefix": pfx})
    # повтор: прошло 60 с < 300 с — цена берётся из кэша (без нового запроса), запуск второй машины
    e.add_response("run_instances", run_resp("i-2"), expected_run_params())
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-2", tags=own_tags())]}]},
                   {"Filters": [{"Name": "instance-id", "Values": ["i-2"]}]})
    done = json.dumps({"phase": "done", "exit_code": 0}).encode()
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    s.add_response("list_objects_v2", {"IsTruncated": False, "Contents": [{"Key": pfx + "_DONE.json"}]},
                   {"Bucket": BUCKET, "Prefix": pfx})
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-2"]})
    e.add_response("describe_instances", {"Reservations": []}, reap_filters())
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res["su2-comp130"]["status"] == "ok" and res["su2-comp130"]["attempts"] == 2
    ends = [x for x in R.Ledger(tmp_path / "ledger.jsonl").records() if x["ev"] == "end"]
    assert [(x["iid"], x["outcome"]) for x in ends] == [("i-1", "lost"), ("i-2", "done")]
    aws.done()


def test_reconcile_closes_dead_instances_from_crashed_session(tmp_path):
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    aws, clk = Aws(), Clock(T0 + 7200)
    aws.st["ec2"].add_response("describe_instances", {"Reservations": []},
                               {"Filters": [{"Name": "instance-id", "Values": ["i-open"]}]})
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    r.reconcile()
    s = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s["open"] == [] and s["committed"] == 0
    # закрыта по верхней границе: (1.00 + 0.01) × 70 мин
    assert s["spent"] == pytest.approx(1.01 + 1.01 * 70 / 60)
    aws.done()


# ----------------------------------------------------------------------------- параллельность и задачи
def test_parallel_limit_and_vcpu_quota(tmp_path):
    def fake_launch(self, job, attempt=1):
        iid = f"i-{job.name}"
        return "launched", R.Active(job, attempt, iid, job.instance_type, "us-east-1a", T0, T0 + 3600)

    jobs = [make_job(name=f"j{i}") for i in range(8)]
    for quota, limit, expect in ((384, 6, 6), (128, 6, 2), (384, 3, 3)):
        r = R.Runner(make_cfg(tmp_path, vcpu_quota=quota, max_parallel=limit), out=lambda m: None)
        r.launch = fake_launch.__get__(r)
        from collections import deque
        pending = deque((j, 1) for j in jobs)
        r._launch_some(pending)
        assert len(r.active) == expect and len(pending) == 8 - expect


def test_load_jobs_validation(tmp_path):
    f = tmp_path / "jobs.json"
    f.write_text(json.dumps({"defaults": {"max_min": 45, "pip": ["numpy"]},
                             "jobs": [{"name": "a", "command": "true"},
                                      {"name": "b", "command": "true", "instance_type": "c7i.8xlarge",
                                       "fallback_types": []}]}))
    jobs = R.load_jobs(f)
    assert [(j.name, j.max_min, j.vcpu, j.types) for j in jobs] == [
        ("a", 45, 64, ["c7a.16xlarge", "c7i.16xlarge"]), ("b", 45, 32, ["c7i.8xlarge"])]
    for bad in ([{"name": "a", "command": "true", "typo": 1}],
                [{"name": "a", "command": "true"}, {"name": "a", "command": "true"}],
                [{"name": "a b", "command": "true"}],
                [{"name": "a", "command": "true", "max_min": 5}],
                [{"name": "a", "command": "true", "instance_type": "x9.huge"}],
                [{"name": "a", "command": "true", "pip": ["numpy; rm -rf /"]}],
                [{"name": "a", "command": "true", "dnf": ["gcc"]}]):
        f.write_text(json.dumps(bad))
        with pytest.raises(R.RunnerError):
            R.load_jobs(f)


def test_job_and_config_guards(tmp_path):
    """Обходы квоты и потолка отвергаются: заниженный/отрицательный vcpu, отрицательный расход и т. п."""
    j = make_job(vcpu=1)                                        # c7a.16xlarge = 64 по таблице
    assert j.vcpu == 64
    with pytest.raises(R.RunnerError):
        make_job(vcpu=-1000)
    with pytest.raises(R.RunnerError):
        make_job(vcpu="64")
    assert make_job(instance_type="m7a.large", fallback_types=[], vcpu=2).vcpu == 2
    f = tmp_path / "j.json"
    for bad in (["oops"], {"defaults": [], "jobs": [{"name": "a", "command": "true"}]},
                [{"name": 5, "command": "true"}], "text"):
        f.write_text(json.dumps(bad))
        with pytest.raises(R.RunnerError):
            R.load_jobs(f)
    f.write_text("{not json")
    with pytest.raises(R.RunnerError):
        R.load_jobs(f)
    for kw in ({"spent_before": -1000.0}, {"grace_min": -30}, {"poll_s": 0.0}, {"retries": -1},
               {"ebs_usd_gb_month": -1.0}, {"ipv4_usd_h": -1.0}, {"price_margin": 0.5}, {"ceiling": 0.0}):
        with pytest.raises(R.RunnerError):
            R.validate_config(make_cfg(tmp_path, **kw))
    R.validate_config(make_cfg(tmp_path))


def test_example_jobs_file_is_valid():
    jobs = R.load_jobs(Path(__file__).resolve().parent / "jobs.example.json")
    assert {j.name for j in jobs} == {"smoke", "explore-249"}


# ----------------------------------------------------------------------------- маскирование и CLI
def test_redactor_and_ledger_guard(tmp_path):
    red = R.Redactor(BUCKET)
    s = red(f"s3://{BUCKET}/x key={FAKE_SECRET} id={FAKE_AKID} other=AKIAPLACEHOLDER00000")
    assert BUCKET not in s and FAKE_SECRET not in s and FAKE_AKID not in s and "AKIAPLACEHOLDER" not in s
    assert "<BUCKET>" in s
    led = R.Ledger(tmp_path / "l.jsonl", red)
    led.append("x", note=f"path s3://{BUCKET}/p", secret=FAKE_SECRET)     # маскируется, не падает
    text = (tmp_path / "l.jsonl").read_text()
    assert BUCKET not in text and FAKE_SECRET not in text
    with pytest.raises(R.RunnerError):                                     # имя переменной ключа — отказ
        led.append("x", note="aws_secret_access_key=...")
    # текст AccessDenied от AWS: номер аккаунта в ARN скрыт
    msg = red("User: arn:aws:iam::123456789012:user/alpha-o-agent is not authorized on "
              "arn:aws:ec2:us-east-1:123456789012:instance/*")
    assert "123456789012" not in msg and msg.count("<аккаунт>") == 2


def test_cli_plan_offline_and_status(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("ALPHA_O_S3", f"s3://{BUCKET}/{PREFIX}")
    jobs = str(Path(__file__).resolve().parent / "jobs.example.json")
    rc = R.main(["plan", jobs, "--assume-price", "1.5", "--show-userdata", "--run-id", "cli-plan",
                 "--ledger", str(tmp_path / "ledger.jsonl"), "--out", str(tmp_path / "runs")])
    out = capsys.readouterr().out
    assert rc == 0 and "shutdown -P +20" in out and '"MarketType": "spot"' in out
    assert BUCKET not in out and FAKE_SECRET not in out and FAKE_AKID not in out
    assert "прогноз пакета 2.41" in out                       # (1.8 + диск) × (20+10 и 40+10 мин)
    assert not (tmp_path / "ledger.jsonl").exists()           # план не пишет в общий журнал
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    assert R.main(["status", "--ledger", str(tmp_path / "ledger.jsonl")]) == 0
    out = capsys.readouterr().out
    assert "израсходовано ≈ 1.01" in out and "по ним 2.34" in out and "остаток 1996.64" in out   # 2.345 → 2.34 (float); потолок по умолчанию 2000 (решение 14.11.2026)


def test_cli_run_refuses_without_env_before_any_aws_call(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("ALPHA_O_S3", raising=False)
    monkeypatch.delenv("ALPHA_O_INSTANCE_PROFILE", raising=False)
    jobs = str(Path(__file__).resolve().parent / "jobs.example.json")
    rc = R.main(["run", jobs, "--dryrun", "--ledger", str(tmp_path / "l.jsonl")])
    err = capsys.readouterr().err
    assert rc == 2 and "ALPHA_O_S3" in err and FAKE_SECRET not in err and FAKE_AKID not in err


def test_cli_run_without_live_is_dryrun(tmp_path, capsys, monkeypatch):
    """`run` без --live: run_instances только с DryRun=True, номер аккаунта в выводе скрыт."""
    monkeypatch.setenv("ALPHA_O_S3", f"s3://{BUCKET}/{PREFIX}")
    monkeypatch.setenv("ALPHA_O_INSTANCE_PROFILE", "gtu-test-profile")
    monkeypatch.delenv("GTU_LIVE", raising=False)
    aws, sts = Aws(), client("sts")
    st_sts = Stubber(sts)
    st_sts.add_response("get_caller_identity", {"UserId": "AIDAEXAMPLE", "Account": "123456789012",
                                                "Arn": "arn:aws:iam::123456789012:user/alpha-o-agent"})
    st_sts.activate()
    aws.st["s3"].add_response("list_objects_v2", {"IsTruncated": False},
                              {"Bucket": BUCKET, "Prefix": f"{PREFIX}/gtu/", "MaxKeys": 1})
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               {"InstanceTypes": ["c7a.16xlarge"], "ProductDescriptions": ["Linux/UNIX"],
                                "StartTime": ANY, "EndTime": ANY})
    add_ami(aws)
    exp = dict(expected_run_params(), DryRun=True, ImageId="ami-0123456789abcdef0")
    exp["TagSpecifications"] = ANY                              # run_id генерируется
    aws.st["ec2"].add_client_error("run_instances", service_error_code="DryRunOperation",
                                   service_message="Request would have succeeded", http_status_code=412,
                                   expected_params=exp)
    monkeypatch.setattr(R, "_clients", lambda region: (aws.ec2, aws.s3, aws.ssm, sts))
    f = tmp_path / "jobs.json"
    f.write_text(json.dumps([{"name": "su2-comp130", "command": "true", "max_min": 60}]))
    rc = R.main(["run", str(f), "--skip-branch-check", "--ledger", str(tmp_path / "ledger.jsonl"),
                 "--out", str(tmp_path / "runs")])
    out = capsys.readouterr().out
    assert rc == 0 and "[умолчание]" in out and "DryRunOperation" in out
    assert "123456789012" not in out and BUCKET not in out
    assert not (tmp_path / "ledger.jsonl").exists()
    aws.done()
    st_sts.assert_no_pending_responses()


def test_default_ceiling_is_owner_decision():
    # Решение владельца 14.11.2026: потолок 2000 долл., отчёт после каждых 500 долл.
    import dataclasses
    f = {x.name: x.default for x in dataclasses.fields(R.Config)}
    assert f["ceiling"] == 2000.0 and f["report_step"] == 500.0
