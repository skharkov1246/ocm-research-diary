"""Тесты runner.py без сети: botocore Stubber и свои заглушки вместо AWS, user-data — на подменённых командах.

  python3 -m pytest -q routes/gtu/platform/aws/test_runner.py

Ключи в тестах поддельные (ставятся в окружение на время теста); любой вызов AWS,
для которого нет заготовленного ответа, Stubber превращает в ошибку — сети нет.
Выдуманные значения помечены словами SYNTHETIC / FIXTURE / EXAMPLE (сканер gate.yml).
"""
import datetime as dt
import io
import json
import os
import re
import shutil
import subprocess
import sys
from collections import deque
from pathlib import Path

import botocore.session
import pytest
from botocore.exceptions import ClientError
from botocore.response import StreamingBody
from botocore.stub import ANY, Stubber

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(HERE))
import runner as R  # noqa: E402

FAKE_AKID = "AKIASYNTHETIC0000000"  # выдуманный ключ: пометка SYNTHETIC по правилу gate.yml
FAKE_SECRET = "fakeSecretValueForTestsOnly/0123456789abcdEF"
BUCKET = "alpha-o-results-123456789012"
PREFIX = "alpha-o"
RUN = "20261114T000000Z-ab12"
T0 = 1_800_000_000.0
UTC = dt.timezone.utc
AMI = "ami-0123456789abcdef0"
# Подсети по умолчанию — выдуманные id (FIXTURE).
SUBNETS = {az: f"subnet-0fixture000000{az[-2:]}" for az in ("us-east-1a", "us-east-1b", "us-east-1c")}
# Цены спота c7a.16xlarge по зонам из живого прогона 14.11.2026 (1a — 0,955; 1b — 0,748; 1c — 0,706).
LIVE_C7A = {"us-east-1a": 0.955, "us-east-1b": 0.748, "us-east-1c": 0.706}


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


def flat(prices):
    """{зона: цена} -> вход spot_resp с одной записью на зону."""
    return {az: [(p, -60)] for az, p in prices.items()}


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


def ledger(tmp_path, ev=None):
    recs = R.Ledger(tmp_path / "ledger.jsonl").records()
    return [x for x in recs if ev is None or x["ev"] == ev]


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


# ----------------------------------------------------------------------------- заготовки run_instances
def expected_run_params(itype="c7a.16xlarge", maxp="1.5600", job="su2-comp130", disk=50, root="/dev/xvda",
                        az="us-east-1a", how="subnet", ami=AMI):
    tags = [{"Key": "Project", "Value": R.DEFAULT_TAG}, {"Key": "RunId", "Value": RUN},
            {"Key": "Job", "Value": job}, {"Key": "Name", "Value": f"gtu-{job}"}]
    p = {
        "ImageId": ami, "InstanceType": itype, "MinCount": 1, "MaxCount": 1,
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
    if how == "subnet":
        p["SubnetId"] = SUBNETS[az]
    elif how == "eni":
        p["NetworkInterfaces"] = [{"DeviceIndex": 0, "SubnetId": SUBNETS[az], "AssociatePublicIpAddress": True,
                                   "DeleteOnTermination": True}]
    elif how == "placement":
        p["Placement"] = {"AvailabilityZone": az}
    return p


def add_ami(aws, image="al2023"):
    aws.st["ssm"].add_response("get_parameter",
                               {"Parameter": {"Name": R.IMAGES[image]["ssm"], "Type": "String",
                                              "Value": AMI, "Version": 1}},
                               {"Name": R.IMAGES[image]["ssm"]})


def add_subnets(aws, azs=tuple(SUBNETS), public=True):
    aws.st["ec2"].add_response("describe_subnets", {"Subnets": [
        {"SubnetId": SUBNETS[az], "AvailabilityZone": az, "MapPublicIpOnLaunch": public, "State": "available"}
        for az in azs]}, {"Filters": [{"Name": "default-for-az", "Values": ["true"]}]})


def run_resp(iid, itype="c7a.16xlarge", az="us-east-1a"):
    return {"Instances": [{"InstanceId": iid, "InstanceType": itype, "Placement": {"AvailabilityZone": az},
                           "State": {"Name": "pending"}}]}


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
    rec = ledger(tmp_path, "budget")[-1]
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
    # run_instances и образ НЕ заготовлены: если исполнитель попробует запустить машину, Stubber упадёт.
    # потолок задан явно: тест проверяет отказ при превышении, а не значение по умолчанию
    r = R.Runner(make_cfg(tmp_path, ceiling=300.0), aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res == {"su2-comp130": {"status": "refused", "attempts": 1}}
    budget = ledger(tmp_path, "budget")
    assert [(b["itype"], b["ok"]) for b in budget] == [("c7a.16xlarge", False), ("c7i.16xlarge", False)]
    assert all(b["total"] > b["ceiling"] for b in budget)
    assert not [x for x in ledger(tmp_path, "launch") if x["run_id"] == RUN]
    aws.done()


def test_ceiling_counts_open_instances_and_parallel_launches(tmp_path):
    """Открытые машины входят в прогноз целиком — параллельный запуск не обходит потолок."""
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    s = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s["open"] == ["i-open"] and s["closed"] == 1
    assert s["spent"] == pytest.approx(1.01) and s["committed"] == pytest.approx(2.345)


def test_second_launch_waits_for_budget_then_starts(tmp_path):
    """Потолок 3 долл.: первая (≈1.83) запущена; вторая с её верхней границей — отказ, но не снята:
    ждёт в очереди и запускается, когда первая закрылась (верхняя граница сменилась фактом)."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-a"), expected_run_params(job="a"))
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", {"us-east-1a": [(1.4, -60)]}),
                   spot_params("c7i.16xlarge"))
    e.add_response("run_instances", run_resp("i-b"), expected_run_params(job="b"))
    r = R.Runner(make_cfg(tmp_path, ceiling=3.0), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    pending = deque([(make_job(name="a"), 1), (make_job(name="b"), 1)])
    assert r._launch_some(pending) is False
    assert list(r.active) == ["i-a"] and "b" not in r.results and pending[0][0].name == "b"
    budget = ledger(tmp_path, "budget")
    assert [(b["job"], b["ok"]) for b in budget] == [("a", True), ("b", False), ("b", False)]
    assert budget[1]["committed"] == pytest.approx(budget[0]["new"], abs=1e-3)
    r._launch_some(pending)                                  # ничего не закрылось — бюджет не пересчитывается
    assert len(ledger(tmp_path, "budget")) == 3
    clk.sleep(60)
    r._close(r.active["i-a"], "done")                        # первая закончилась за 60 с
    r._launch_some(pending)
    assert list(r.active) == ["i-b"] and not pending
    assert [(b["job"], b["ok"]) for b in ledger(tmp_path, "budget")][-1] == ("b", True)
    aws.done()


# ----------------------------------------------------------------------------- зоны (пункт А)
def test_cheapest_zone_first_subnet_and_ledger(tmp_path):
    """Зоны по возрастанию цены: машина просится в 1c (0,706), подсеть 1c; бюджет и MaxPrice — по максимуму."""
    aws, clk, out = Aws(), Clock(), []
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-c", az="us-east-1c"),
                   expected_run_params(az="us-east-1c", maxp="1.1460"))   # ceil(0.955 × 1.2) — максимум, не 1c
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=out.append)
    status, act = r.launch(make_job())
    assert status == "launched" and act.az == "us-east-1c" and act.az_rank == 0 and act.az_cheapest == "us-east-1c"
    ud = r.last_params["UserData"]
    assert "GTU_AZ=us-east-1c GTU_AZ_CHEAPEST=us-east-1c" in ud
    rec = ledger(tmp_path, "launch")[0]
    assert (rec["az"], rec["az_req"], rec["az_rank"], rec["spot_usd_h"]) == ("us-east-1c", "us-east-1c", 0, 0.706)
    assert rec["zones"] == [["us-east-1c", 0.706], ["us-east-1b", 0.748], ["us-east-1a", 0.955]]
    assert rec["max_usd_h"] == pytest.approx(1.146) and rec["net"] == "subnet"
    assert not any(s in json.dumps(rec) for s in SUBNETS.values())          # id подсети в журнал не попадает
    b = ledger(tmp_path, "budget")[0]
    assert b["max_usd_h"] == pytest.approx(1.146) and b["spot_usd_h"] == pytest.approx(0.955)
    assert any("us-east-1c 0.7060 → us-east-1b 0.7480 → us-east-1a 0.9550" in m for m in out)
    aws.done()


def test_zone_fallback_then_type_fallback(tmp_path):
    """1c — нет ёмкости, 1b — SpotMaxPriceTooLow, 1a — Unsupported → запасной c7i в его самой дешёвой зоне."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    for az, code in (("us-east-1c", "InsufficientInstanceCapacity"), ("us-east-1b", "SpotMaxPriceTooLow"),
                     ("us-east-1a", "Unsupported")):
        e.add_client_error("run_instances", service_error_code=code, http_status_code=500,
                           expected_params=expected_run_params(az=az, maxp="1.1460"))
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", flat({"us-east-1b": 0.80, "us-east-1c": 0.78})),
                   spot_params("c7i.16xlarge"))
    e.add_response("run_instances", run_resp("i-i", "c7i.16xlarge", "us-east-1c"),
                   expected_run_params(itype="c7i.16xlarge", az="us-east-1c", maxp="0.9600"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and (act.itype, act.az, act.az_rank) == ("c7i.16xlarge", "us-east-1c", 0)
    errs = [(x["itype"], x["az"], x["code"]) for x in ledger(tmp_path, "launch_error")]
    assert errs == [("c7a.16xlarge", "us-east-1c", "InsufficientInstanceCapacity"),
                    ("c7a.16xlarge", "us-east-1b", "SpotMaxPriceTooLow"),
                    ("c7a.16xlarge", "us-east-1a", "Unsupported")]
    # зоны со сбоем ёмкости — в конце очереди для следующих задач (на avoid_s)
    assert [z for z, _ in r.zone_order("c7a.16xlarge", r.price("c7a.16xlarge"))][-1] in ("us-east-1a", "us-east-1b", "us-east-1c")
    assert all(r._avoided("c7a.16xlarge", z) for z in LIVE_C7A)
    aws.done()


def test_fallback_type_on_capacity_error(tmp_path):
    """Старый сценарий: обе зоны c7a без ёмкости — запасной тип со своей проверкой бюджета."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    for az in ("us-east-1a", "us-east-1b"):
        e.add_client_error("run_instances", service_error_code="InsufficientInstanceCapacity",
                           service_message="no capacity", http_status_code=500,
                           expected_params=expected_run_params(az=az))
    e.add_response("describe_spot_price_history", spot_resp("c7i.16xlarge", {"us-east-1c": [(1.0, -60)]}),
                   spot_params("c7i.16xlarge"))
    e.add_response("run_instances", run_resp("i-0bbb", "c7i.16xlarge", "us-east-1c"),
                   expected_run_params(itype="c7i.16xlarge", maxp="1.2000", az="us-east-1c"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and act.itype == "c7i.16xlarge"
    evs = [(x["ev"], x.get("itype"), x.get("code")) for x in ledger(tmp_path)]
    assert ("launch_error", "c7a.16xlarge", "InsufficientInstanceCapacity") in evs
    assert ("launch", "c7i.16xlarge", None) in evs
    aws.done()


def test_subnets_denied_falls_back_to_placement(tmp_path):
    """Нет прав на DescribeSubnets — зона задаётся Placement.AvailabilityZone (подсеть выберет EC2)."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    add_ami(aws)
    e.add_client_error("describe_subnets", service_error_code="UnauthorizedOperation", http_status_code=403,
                       expected_params={"Filters": [{"Name": "default-for-az", "Values": ["true"]}]})
    e.add_response("run_instances", run_resp("i-p", az="us-east-1c"),
                   expected_run_params(az="us-east-1c", maxp="1.1460", how="placement"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    assert r.launch(make_job())[0] == "launched"
    assert "SubnetId" not in r.last_params
    assert ledger(tmp_path, "subnets")[0]["code"] == "UnauthorizedOperation"
    assert ledger(tmp_path, "launch")[0]["net"] == "placement"
    aws.done()


def test_private_default_subnet_gets_explicit_public_ip(tmp_path):
    """Подсеть без MapPublicIpOnLaunch: сетевой интерфейс с публичным IPv4 — машине нужен интернет."""
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws, public=False)
    e.add_response("run_instances", run_resp("i-n", az="us-east-1c"),
                   expected_run_params(az="us-east-1c", maxp="1.1460", how="eni"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    assert r.launch(make_job())[0] == "launched"
    aws.done()


def test_zone_without_default_subnet_is_skipped(tmp_path):
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws, azs=("us-east-1a", "us-east-1b"))
    e.add_response("run_instances", run_resp("i-s", az="us-east-1b"),
                   expected_run_params(az="us-east-1b", maxp="1.1460"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and act.az == "us-east-1b" and act.az_rank == 1 and act.az_cheapest == "us-east-1c"
    aws.done()


# ----------------------------------------------------------------------------- образ (пункт Б)
def test_ubuntu_image_via_describe_images_when_ssm_denied(tmp_path):
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    aws.st["ssm"].add_client_error("get_parameter", service_error_code="AccessDeniedException", http_status_code=400,
                                   expected_params={"Name": R.IMAGES["ubuntu2404"]["ssm"]})
    img = R.IMAGES["ubuntu2404"]
    e.add_response("describe_images", {"Images": [
        {"ImageId": "ami-0fixture0000000old", "CreationDate": "2026-10-01T00:00:00.000Z", "Name": "noble-old-FIXTURE"},
        {"ImageId": "ami-0fixture0000000new", "CreationDate": "2026-11-01T00:00:00.000Z", "Name": "noble-new-FIXTURE"}]},
        {"Owners": [img["owner"]], "Filters": [{"Name": "name", "Values": [img["name"]]},
                                               {"Name": "state", "Values": ["available"]},
                                               {"Name": "architecture", "Values": ["x86_64"]}]})
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-u", az="us-east-1c"),
                   expected_run_params(az="us-east-1c", maxp="1.1460", root="/dev/sda1", ami="ami-0fixture0000000new"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    assert r.launch(make_job(image="ubuntu2404", apt=["gmsh"]))[0] == "launched"
    aws.done()


def test_unavailable_image_fails_only_its_job(tmp_path):
    aws, clk = Aws(), Clock()
    e = aws.st["ec2"]
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", flat(LIVE_C7A)), spot_params("c7a.16xlarge"))
    aws.st["ssm"].add_client_error("get_parameter", service_error_code="AccessDeniedException", http_status_code=400,
                                   expected_params={"Name": R.IMAGES["ubuntu2404"]["ssm"]})
    e.add_client_error("describe_images", service_error_code="UnauthorizedOperation", http_status_code=403)
    add_ami(aws, "al2023")
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-ok", az="us-east-1c"),
                   expected_run_params(job="b", az="us-east-1c", maxp="1.1460"))
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    pending = deque([(make_job(name="a", apt=["gmsh"]), 1), (make_job(name="b"), 1)])
    r._launch_some(pending)
    assert r.results["a"]["status"] == "image_unavailable"
    assert "al2023" in r.results["a"]["detail"] and "container" in r.results["a"]["detail"]
    assert list(r.active) == ["i-ok"]
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
    assert "@@" not in ud                                       # все подстановки сделаны


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


def test_userdata_python_mode_lines(tmp_path):
    job = make_job(python="3.12", pip=["cantera", "numpy"], requirements="routes/gtu/req.txt")
    ud = R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket=BUCKET,
                          key_prefix="p/gtu/r/j", az="us-east-1c", az_cheapest="us-east-1c", itype="c7a.16xlarge")
    (tmp_path / "ud.sh").write_text(ud)
    assert subprocess.run(["bash", "-n", str(tmp_path / "ud.sh")]).returncode == 0
    for s in ('timeout 20m dnf -q install -y git python3-pip || exit 11',       # al2023 по умолчанию
              'python3 -m venv "$G/uvboot" && timeout 10m "$G/uvboot/bin/pip" install -q uv',
              'releases/latest/download/uv-$(uname -m)-unknown-linux-gnu.tar.gz',
              'timeout 15m "$UV" python install 3.12 || exit 19',
              'timeout 10m "$UV" venv --seed -q --python 3.12 "$G/venv" || exit 12',
              'timeout 20m "$UV" pip install -q --python "$G/venv/bin/python" boto3 || exit 13',
              'timeout 40m "$UV" pip install -q --python "$G/venv/bin/python" cantera numpy || exit 17',
              'timeout 40m "$UV" pip install -q --python "$G/venv/bin/python" -r routes/gtu/req.txt || exit 23',
              'export PATH="$G/venv/bin:$PATH" OMP_NUM_THREADS="$NPROC"', "NPROC=$(nproc)", "GTU_ITYPE=c7a.16xlarge"):
        assert s in ud, s
    assert ud.index("phase clone") < ud.index("phase requirements") < ud.index("phase command")
    assert "docker" not in ud


def test_userdata_container_mode_lines(tmp_path):
    job = make_job(container="ubuntu:24.04", container_setup=["apt-get update -q", "apt-get install -y -q calculix-ccx"])
    ud = R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket=BUCKET,
                          key_prefix="p/gtu/r/j")
    (tmp_path / "ud.sh").write_text(ud)
    assert subprocess.run(["bash", "-n", str(tmp_path / "ud.sh")]).returncode == 0
    for s in ("timeout 20m dnf -q install -y git python3-pip docker || exit 11",
              "timeout 5m systemctl start docker || exit 21",
              "timeout 30m docker pull -q ubuntu:24.04 || exit 22",
              "docker run --rm --name gtu-job --user 0:0 --entrypoint /bin/bash",
              '--cpus "$NPROC" --shm-size "$((MEM_KB / 2))k"',
              '-v "$G/repo:/gtu/repo" -v "$OUT:/gtu/out"', "-e GTU_OUT=/gtu/out", "-e NPROC",
              "ubuntu:24.04 /gtu/ctr/entry.sh", "apt-get install -y -q calculix-ccx",
              'cat >"$G/ctr/cmd.sh" <<\'GTU_CMD_EOF\''):
        assert s in ud, s
    fin = ud[ud.index("finish() {"):ud.index("trap finish EXIT")]
    assert fin.index("docker kill gtu-job") < fin.index('up "$OUT"')     # контейнер снят до выгрузки
    assert "uvboot" not in ud and "venv/bin/pip\" install -q --upgrade pip boto3" in ud   # boto3 — снаружи


def test_secrets_refused():
    with pytest.raises(R.RunnerError):
        make_job(command=f"export AWS_ACCESS_KEY_ID={FAKE_AKID}; run")
    with pytest.raises(R.RunnerError):
        make_job(command=f"echo {FAKE_SECRET}")
    with pytest.raises(R.RunnerError):
        make_job(command="echo AKIAEXAMPLE000000000")
    with pytest.raises(R.RunnerError):
        make_job(container="ubuntu:24.04", container_setup=[f"echo {FAKE_SECRET}"])
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
    "dnf": 'echo "dnf $*" >> "$CALLS"\ncase " $* " in *" failpkg "*) exit 1;; esac\n',
    "systemctl": 'echo "systemctl $*" >> "$CALLS"\n',
    "git": ('echo "git $*" >> "$CALLS"\nfor last; do :; done\n'
            'if [ "$1" = clone ]; then mkdir -p "$last"; [ -n "${FAKE_REQ:-}" ] && echo numpy > "$last/$FAKE_REQ"; fi\n'
            'exit 0\n'),
    # IMDS: токен и зона — выдуманные (FIXTURE); прочие адреса недоступны
    "curl": ('case "$*" in\n  *api/token*) echo TOKEN-FIXTURE;;\n  *placement/availability-zone*) echo us-east-1c;;\n'
             '  *instance-type*) echo c7a.16xlarge;;\n  *) exit 22;;\nesac\n'),
    "python3": ('if [ "$1" = -m ] && [ "$2" = venv ]; then\n'
                '  mkdir -p "$3/bin"\n'
                "  printf '#!/bin/bash\\nexec env PYTHONPATH=\"$FAKEMODS\" \"$REALPY\" \"$@\"\\n' > \"$3/bin/python\"\n"
                "  printf '#!/bin/bash\\necho \"pip $*\" >> \"$CALLS\"\\n' > \"$3/bin/pip\"\n"
                '  chmod +x "$3/bin/python" "$3/bin/pip"; echo venv >> "$CALLS"; exit 0\n'
                'fi\nexec "$REALPY" "$@"\n'),
    "uv": ('echo "uv $*" >> "$CALLS"\n'
           'if [ "$1" = venv ]; then\n  for last; do :; done\n  mkdir -p "$last/bin"\n'
           "  printf '#!/bin/bash\\nexec env PYTHONPATH=\"$FAKEMODS\" \"$REALPY\" \"$@\"\\n' > \"$last/bin/python\"\n"
           "  printf '#!/bin/bash\\necho \"pip $*\" >> \"$CALLS\"\\n' > \"$last/bin/pip\"\n"
           '  chmod +x "$last/bin/python" "$last/bin/pip"\nfi\n'
           '[ "$1" = --version ] && echo "uv 0.0.0-FIXTURE"\nexit 0\n'),
    # docker: тома -v — символьные ссылки, -e — окружение, команда — bash на хосте (точка входа /bin/bash)
    "docker": ('echo "docker $*" >> "$CALLS"\n'
               'case "$1" in\n'
               '  pull) [ -n "${FAIL_PULL:-}" ] && exit 1; exit 0;;\n'
               '  run) shift; envs=(); wd=/\n'
               '    while [ $# -gt 0 ]; do\n'
               '      case "$1" in\n'
               '        -v) src=${2%%:*}; rest=${2#*:}; dst=${rest%%:*}; mkdir -p "$(dirname "$dst")"; ln -sfn "$src" "$dst"; shift 2;;\n'
               '        -e) case "$2" in *=*) envs+=("$2");; *) envs+=("$2=${!2}");; esac; shift 2;;\n'
               '        -w) wd=$2; shift 2;;\n'
               '        --entrypoint|--name|--cpus|--shm-size|--user) shift 2;;\n'
               '        --rm|--init) shift;;\n'
               '        *) break;;\n'
               '      esac\n'
               '    done\n'
               '    shift; cd "$wd" || exit 125; exec env "${envs[@]}" bash "$@";;\n'
               'esac\nexit 0\n'),
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


def run_box(tmp_path, job, env_extra=None, skip=(), **ud_kw):
    """Выполнить user-data локально: системные команды подменены, boto3 — поддельный."""
    shims, mods, s3dir, calls = tmp_path / "shims", tmp_path / "mods" / "boto3", tmp_path / "s3", tmp_path / "calls"
    shims.mkdir(parents=True)
    mods.mkdir(parents=True)
    for name, text in SHIMS.items():
        if name in skip:
            continue
        (shims / name).write_text("#!/bin/bash\n" + text)
        (shims / name).chmod(0o755)
    (mods / "__init__.py").write_text(FAKE_BOTO3)
    kp = f"pfx/gtu/rid1/{job.name}"
    ud = R.build_userdata(job, run_id="rid1", attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket="bkt",
                          key_prefix=kp, _root=str(tmp_path / "box"), _log=str(tmp_path / "box.log"),
                          _ctr_root=str(tmp_path / "ctrroot"), **ud_kw)
    (tmp_path / "ud.sh").write_text(ud)
    path = os.environ["PATH"]
    if "uv" in skip:                                            # настоящий uv в PATH не должен подхватиться
        path = ":".join(p for p in path.split(":") if not (Path(p) / "uv").exists())
    env = {"PATH": f"{shims}:{path}", "CALLS": str(calls), "S3DIR": str(s3dir),
           "FAKEMODS": str(tmp_path / "mods"), "REALPY": sys.executable, "HOME": str(tmp_path),
           "REALTIMEOUT": shutil.which("timeout")}
    env.update(env_extra or {})
    subprocess.run(["bash", str(tmp_path / "ud.sh")], env=env, timeout=120, check=False)
    return calls.read_text().splitlines(), s3dir / "bkt" / kp


def _done(s3):
    d = json.loads((s3 / "_DONE.json").read_text())
    return d["phase"], d["exit_code"]


def test_userdata_runs_locally_success(tmp_path):
    job = make_job(name="box-ok", max_min=30, pip=["cantera"],
                   command='echo hello > "$GTU_OUT/result.txt"\nmkdir -p "$GTU_OUT/sub"\n'
                           'echo 42 > "$GTU_OUT/sub/x.txt"\npwd > "$GTU_OUT/pwd.txt"\necho "$NPROC" > "$GTU_OUT/nproc.txt"')
    calls, s3 = run_box(tmp_path, job, az="us-east-1c", az_cheapest="us-east-1c")
    assert calls[0] == "shutdown -P +30"                       # сторож первым
    assert calls[-1] == "poweroff"                              # гасится сама
    assert calls[-2] == "upload pfx/gtu/rid1/box-ok/_DONE.json" # метка — последней
    assert any(c.startswith("git clone -q --depth 1 -b main https://github.com/") for c in calls)
    assert "dnf -q install -y git python3-pip" in calls         # al2023 — образ по умолчанию
    assert "pip install -q cantera" in calls
    assert (s3 / "result.txt").read_text() == "hello\n"
    assert (s3 / "sub" / "x.txt").read_text() == "42\n"
    assert (s3 / "pwd.txt").read_text().strip() == str(tmp_path / "box" / "repo")
    assert (s3 / "nproc.txt").read_text().strip() == subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip()
    assert (s3 / "_job.log").exists() and (s3 / "_stdout.log").exists()
    done = json.loads((s3 / "_DONE.json").read_text())
    assert done["phase"] == "done" and done["exit_code"] == 0 and done["run_id"] == "rid1"
    assert done["az"] == "us-east-1c" and done["az_planned"] == "us-east-1c" and done["itype"] == "c7a.16xlarge"
    assert done["t_start"] <= done["t_cmd"] <= done["t_end"]
    log = (s3 / "_job.log").read_text()
    for ph in ("base", "venv", "packages", "clone", "command"):  # время каждой фазы — в журнале машины
        assert re.search(rf"\[gtu-time\] {ph} \d+s total=\d+s", log), ph
    # выгрузка укладывается до сторожа: результаты — не дольше (срок − прошло − 120 с), метка — 60 с
    ups = [c.split()[1] for c in calls if c.startswith("timeout ") and "upload.py" in c]
    assert len(ups) == 2 and ups[1] == "60s"
    assert 30 * 60 - 120 - 60 <= int(ups[0].rstrip("s")) <= 30 * 60 - 120


def test_upload_window_fits_inside_watchdog():
    """Худший случай: команда снята по timeout -k 60 в последний момент — метка всё равно до сторожа."""
    for max_min in (20, 60, 720):
        for kw in ({}, {"python": "3.12"}, {"container": "ubuntu:24.04"}):
            job = make_job(max_min=max_min, **kw)
            ud = R.build_userdata(job, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main",
                                  bucket=BUCKET, key_prefix="p/gtu/r/j")
            reserve = 8 * 60
            cmd_end = max_min * 60 - reserve + 60                   # DEADLINE + -k 60 (от T0)
            left = min(1800, max(30, max_min * 60 - cmd_end - 120))  # как в finish()
            assert cmd_end + left + 60 <= max_min * 60 - 60         # результаты + метка, минута до сторожа
            assert 'up "$OUT" "$LEFT"' in ud and 'up "$G/_DONE.json" 60' in ud
            assert f"WD_END=$((T0 + {max_min} * 60))" in ud


def test_userdata_runs_locally_failures(tmp_path):
    calls, s3 = run_box(tmp_path / "a", make_job(name="box-fail", max_min=30, command="exit 3"))
    assert _done(s3) == ("done", 3)
    assert calls[0] == "shutdown -P +30" and calls[-1] == "poweroff"
    # сбой пакетов задачи: журнал и метка всё равно выгружены, клонирования нет
    calls, s3 = run_box(tmp_path / "b", make_job(name="box-pkg", max_min=30, apt=["failpkg"]))
    assert _done(s3) == ("packages", 16)
    assert not any(c.startswith("git clone") for c in calls)
    assert (s3 / "_job.log").exists() and calls[-1] == "poweroff"


def test_userdata_python_mode_runs_locally(tmp_path):
    job = make_job(name="box-py", max_min=30, python="3.12", pip=["cantera"], requirements="req.txt",
                   command='command -v python > "$GTU_OUT/which.txt"\necho "$NPROC" > "$GTU_OUT/nproc.txt"')
    calls, s3 = run_box(tmp_path, job, env_extra={"FAKE_REQ": "req.txt"})
    venv = tmp_path / "box" / "venv"
    assert _done(s3) == ("done", 0) and calls[-1] == "poweroff"
    assert "uv python install 3.12" in calls
    assert f"uv venv --seed -q --python 3.12 {venv}" in calls
    assert f"uv pip install -q --python {venv}/bin/python boto3" in calls
    assert f"uv pip install -q --python {venv}/bin/python cantera" in calls
    assert f"uv pip install -q --python {venv}/bin/python -r req.txt" in calls
    assert calls.index(f"uv pip install -q --python {venv}/bin/python -r req.txt") > \
        next(i for i, c in enumerate(calls) if c.startswith("git clone"))        # requirements — после клона
    assert (s3 / "which.txt").read_text().strip() == f"{venv}/bin/python"        # PATH — на venv
    assert int((s3 / "nproc.txt").read_text()) >= 1
    log = (s3 / "_job.log").read_text()
    assert "uv 0.0.0-FIXTURE" in log and re.search(r"\[gtu-time\] venv \d+s", log)
    assert re.search(r"\[gtu-time\] requirements \d+s", log)


def test_userdata_python_mode_without_uv_fails_cleanly(tmp_path):
    """uv не поставился ни из PyPI, ни из выпуска GitHub: фаза venv, код 18, журнал выгружен, машина гаснет."""
    job = make_job(name="box-nouv", max_min=30, python="3.12")
    calls, s3 = run_box(tmp_path, job, skip=("uv",))
    assert "pip install -q uv" in calls                          # сначала PyPI в отдельный venv
    assert calls[-1] == "poweroff"
    # boto3 в venv задачи не поставлен — выгрузка через aws cli или «UPLOAD-UNAVAILABLE»: метки в S3 может не быть
    log = (tmp_path / "box.log").read_text()
    assert "finish phase=venv rc=18" in log and "GitHub release" in log


def _ctr_job(**kw):
    base = dict(name="box-ctr", max_min=30, container="ubuntu:24.04", container_setup=["echo setup-ok"],
                command='echo "inside nproc=$NPROC" > "$GTU_OUT/inside.txt"\npwd > "$GTU_OUT/pwd.txt"')
    base.update(kw)
    return make_job(**base)


def test_userdata_container_runs_locally(tmp_path):
    calls, s3 = run_box(tmp_path, _ctr_job())
    assert _done(s3) == ("done", 0)
    assert calls[0] == "shutdown -P +30" and calls[-1] == "poweroff"
    assert calls[-2] == "upload pfx/gtu/rid1/box-ctr/_DONE.json"
    assert "dnf -q install -y git python3-pip docker" in calls and "systemctl start docker" in calls
    assert "docker pull -q ubuntu:24.04" in calls
    run = next(c for c in calls if c.startswith("docker run"))
    assert "--cpus" in run and "--shm-size" in run and "-e NPROC" in run
    assert "docker kill gtu-job" in calls                         # снаружи, в finish
    nproc = subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip()
    assert (s3 / "inside.txt").read_text().strip() == f"inside nproc={nproc}"
    assert (s3 / "pwd.txt").read_text().strip() == str(tmp_path / "ctrroot" / "repo")
    assert (s3 / "_container_setup.log").read_text().strip() == "setup-ok"
    assert re.search(r"\[gtu-time\] container_setup \d+s", (s3 / "_job.log").read_text())


@pytest.mark.parametrize("case,kw,env,expect", [
    ("setup-fails", {"container_setup": ["false"]}, {}, ("container", 20)),
    ("command-fails", {"command": "exit 7"}, {}, ("done", 7)),
    ("pull-fails", {}, {"FAIL_PULL": "1"}, ("pull", 22)),
])
def test_container_failure_still_uploads_and_powers_off(tmp_path, case, kw, env, expect):
    calls, s3 = run_box(tmp_path, _ctr_job(**kw), env_extra=env)
    assert _done(s3) == expect
    assert calls[-1] == "poweroff" and calls[-2].endswith("/_DONE.json")
    assert (s3 / "_job.log").exists()
    if case == "setup-fails":
        assert not (s3 / "inside.txt").exists()
    if case == "pull-fails":
        assert not any(c.startswith("git clone") for c in calls) and not any(c.startswith("docker run") for c in calls)


# ----------------------------------------------------------------------------- run_instances
def test_run_instances_params(tmp_path):
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    aws.st["ec2"].add_response("run_instances", run_resp("i-0aaa"), expected_run_params())
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, out=lambda m: None)
    status, act = r.launch(make_job())
    assert status == "launched" and act.iid == "i-0aaa" and act.az == "us-east-1a"
    assert act.deadline == T0 + 70 * 60
    p = r.last_params
    assert p["InstanceMarketOptions"]["MarketType"] == "spot"
    assert p["MetadataOptions"]["HttpTokens"] == "required"
    assert p["InstanceInitiatedShutdownBehavior"] == "terminate"
    assert p["UserData"].splitlines()[1] == "shutdown -P +60"
    assert len(p["ClientToken"]) <= 64 and "KeyName" not in p
    rec = ledger(tmp_path, "launch")[0]
    assert rec["iid"] == "i-0aaa" and rec["spot_usd_h"] == 1.1 and rec["max_usd_h"] == pytest.approx(1.56)
    assert rec["allot_min"] == 70 and rec["az_rank"] == 0 and rec["image"] == "al2023" and rec["mode"] == "system"
    aws.done()


def test_dryrun_does_not_launch(tmp_path):
    aws, clk = Aws(), Clock()
    aws.st["ec2"].add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A),
                               spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
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
    assert any(x["ev"] == "dryrun" and x["ok"] and x["az"] == "us-east-1a" for x in side)
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


def test_attempt_prefixes_do_not_collide_with_jobs(tmp_path):
    r = R.Runner(make_cfg(tmp_path), out=lambda m: None)
    assert r.key_prefix("a") == f"{PREFIX}/gtu/{RUN}/a"
    assert r.key_prefix("a", attempt=2) == f"{PREFIX}/gtu/{RUN}/_attempt2/a"
    with pytest.raises(R.RunnerError):                          # задача не может называться «_attempt2»
        make_job(name="_attempt2")


# ----------------------------------------------------------------------------- полный цикл
def test_full_cycle_launch_poll_collect_reap_no_secrets(tmp_path):
    aws, clk, out = Aws(), Clock(), []
    e, s = aws.st["ec2"], aws.st["s3"]
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-0ccc"), expected_run_params())
    filt = {"Filters": [{"Name": "instance-id", "Values": ["i-0ccc"]}]}
    # опрос 1: машина работает, метки нет
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-0ccc", tags=own_tags())]}]}, filt)
    s.add_client_error("get_object", service_error_code="NoSuchKey", http_status_code=404,
                       expected_params={"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    # опрос 2: метка есть -> сбор -> гашение
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-0ccc", tags=own_tags())]}]}, filt)
    done = json.dumps({"run_id": RUN, "job": "su2-comp130", "attempt": 1, "phase": "done", "exit_code": 0,
                       "az": "us-east-1a"}).encode()
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
    j = res["su2-comp130"]
    assert j["status"] == "ok" and j["files"] == 2
    assert (j["az"], j["az_rank"], j["cheapest_zone"], j["az_machine"]) == ("us-east-1a", 0, True, "us-east-1a")
    assert (tmp_path / "runs" / RUN / "su2-comp130" / "a.txt").read_text() == "hi\n"
    s_ = R.Ledger(tmp_path / "ledger.jsonl").summary()
    assert s_["open"] == [] and s_["closed"] == 1
    assert s_["spent"] == pytest.approx((1.1 + 50 * 0.08 / 730 + 0.005) * 120 / 3600, abs=1e-4)   # 2 опроса по 60 с
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
    add_subnets(aws)
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
    """Машина погашена без метки (прерывание спота) -> одна повторная попытка в свой префикс и в другую зону."""
    aws, clk = Aws(), Clock()
    e, s = aws.st["ec2"], aws.st["s3"]
    pfx = f"{PREFIX}/gtu/{RUN}/su2-comp130/"
    pfx2 = f"{PREFIX}/gtu/{RUN}/_attempt2/su2-comp130/"
    e.add_response("describe_spot_price_history", spot_resp("c7a.16xlarge", PRICES_C7A), spot_params("c7a.16xlarge"))
    add_ami(aws)
    add_subnets(aws)
    e.add_response("run_instances", run_resp("i-1", az="us-east-1a"), expected_run_params())
    e.add_response("describe_instances", {"Reservations": [{"Instances": [
        inst("i-1", "terminated", own_tags(), az="us-east-1a", reason="Server.SpotInstanceTermination")]}]},
        {"Filters": [{"Name": "instance-id", "Values": ["i-1"]}]})
    s.add_client_error("get_object", service_error_code="NoSuchKey", http_status_code=404,
                       expected_params={"Bucket": BUCKET, "Key": pfx + "_DONE.json"})
    s.add_response("list_objects_v2", {"IsTruncated": False}, {"Bucket": BUCKET, "Prefix": pfx})
    # повтор: цена из кэша (60 с < 300 с); зона 1a после прерывания — в конце, просимся в 1b
    e.add_response("run_instances", run_resp("i-2", az="us-east-1b"), expected_run_params(az="us-east-1b"))
    e.add_response("describe_instances", {"Reservations": [{"Instances": [inst("i-2", tags=own_tags())]}]},
                   {"Filters": [{"Name": "instance-id", "Values": ["i-2"]}]})
    done = json.dumps({"phase": "done", "exit_code": 0}).encode()
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx2 + "_DONE.json"})
    s.add_response("list_objects_v2", {"IsTruncated": False, "Contents": [{"Key": pfx2 + "_DONE.json"}]},
                   {"Bucket": BUCKET, "Prefix": pfx2})
    s.add_response("get_object", {"Body": body(done)}, {"Bucket": BUCKET, "Key": pfx2 + "_DONE.json"})
    e.add_response("terminate_instances", {"TerminatingInstances": []}, {"InstanceIds": ["i-2"]})
    e.add_response("describe_instances", {"Reservations": []}, reap_filters())
    r = R.Runner(make_cfg(tmp_path), aws.ec2, aws.s3, aws.ssm, clock=clk, sleep=clk.sleep, out=lambda m: None)
    res = r.run([make_job()])
    assert res["su2-comp130"]["status"] == "ok" and res["su2-comp130"]["attempts"] == 2
    assert res["su2-comp130"]["az"] == "us-east-1b" and res["su2-comp130"]["az_rank"] == 1
    ends = [x for x in ledger(tmp_path, "end")]
    assert [(x["iid"], x["outcome"]) for x in ends] == [("i-1", "lost"), ("i-2", "done")]
    assert "GTU_ATTEMPT=2" in r.last_params["UserData"] and f"S3_KEY={pfx2.rstrip('/')}" in r.last_params["UserData"]
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


# ----------------------------------------------------------------------------- очередь (пункт Г): своя заглушка AWS
def cerr(code, op="Op"):
    return ClientError({"Error": {"Code": code, "Message": code}}, op)


class FakeCloud:
    """Своя заглушка EC2/S3/SSM: машина «работает» run_s секунд и выгружает результат и метку, как user-data;
    прерывание спота — по таблице interrupt {(задача, попытка): через сколько секунд}."""

    def __init__(self, clock, run_s=300, prices=None, fail=(), interrupt=None, irq_reason=R.SPOT_REASON):
        self.clock, self.run_s, self.irq_reason = clock, run_s, irq_reason
        self.prices = prices or {"c7a.16xlarge": LIVE_C7A, "c7i.16xlarge": {"us-east-1b": 0.80, "us-east-1c": 0.78}}
        self.fail, self.interrupt = list(fail), dict(interrupt or {})
        self.inst, self.objects, self.runs, self.max_live, self.max_vcpu = {}, {}, [], 0, 0
        self.ec2, self.s3, self.ssm = _FakeEc2(self), _FakeS3(self), _FakeSsm()

    def live(self):
        return [i for i in self.inst.values() if i["state"] in ("pending", "running")]

    def tick(self):
        now = self.clock()
        for i in self.inst.values():
            if i["state"] != "running" or i["uploaded"]:
                continue
            irq = self.interrupt.get((i["job"], i["attempt"]))
            if irq is not None and now - i["t"] >= irq:
                i["state"], i["reason"], i["uploaded"] = "terminated", self.irq_reason, True
                self.objects[i["key"] + "/partial.txt"] = b"partial"
                self.objects[i["key"] + "/_DONE.json"] = json.dumps(
                    {"attempt": i["attempt"], "phase": "command", "exit_code": 143}).encode()
            elif now - i["t"] >= self.run_s:
                i["uploaded"] = True
                self.objects[i["key"] + "/result.txt"] = f"{i['job']} a{i['attempt']}".encode()
                self.objects[i["key"] + "/_DONE.json"] = json.dumps(
                    {"attempt": i["attempt"], "phase": "done", "exit_code": 0, "az": i["az"]}).encode()


class _FakeEc2:
    def __init__(self, c):
        self.c = c

    def describe_spot_price_history(self, InstanceTypes, **kw):
        t = dt.datetime.fromtimestamp(self.c.clock() - 60, tz=UTC)
        return {"SpotPriceHistory": [{"AvailabilityZone": az, "InstanceType": InstanceTypes[0], "SpotPrice": f"{p:.4f}",
                                      "Timestamp": t} for az, p in self.c.prices.get(InstanceTypes[0], {}).items()]}

    def describe_subnets(self, **kw):
        return {"Subnets": [{"SubnetId": s, "AvailabilityZone": az, "MapPublicIpOnLaunch": True, "State": "available"}
                            for az, s in SUBNETS.items()]}

    def run_instances(self, **p):
        self.c.runs.append(p)
        if self.c.fail:
            code = self.c.fail.pop(0)
            if code:
                raise cerr(code, "RunInstances")
        az = next(a for a, s in SUBNETS.items() if s == p["SubnetId"])
        tags = {t["Key"]: t["Value"] for t in p["TagSpecifications"][0]["Tags"]}
        key = re.search(r"S3_KEY=(\S+)", p["UserData"]).group(1).strip("'")
        attempt = int(re.search(r"GTU_ATTEMPT=(\d+)", p["UserData"]).group(1))
        iid = f"i-{len(self.c.inst) + 1:04d}"
        self.c.inst[iid] = dict(state="running", t=self.c.clock(), job=tags["Job"], attempt=attempt, az=az, key=key,
                                tags=tags, reason="", uploaded=False, vcpu=R.VCPU[p["InstanceType"]],
                                token=p["ClientToken"])
        live = self.c.live()
        self.c.max_live = max(self.c.max_live, len(live))
        self.c.max_vcpu = max(self.c.max_vcpu, sum(i["vcpu"] for i in live))
        return {"Instances": [{"InstanceId": iid, "Placement": {"AvailabilityZone": az}}]}

    def describe_instances(self, Filters, **kw):
        self.c.tick()
        f = {x["Name"]: x["Values"] for x in Filters}
        out = []
        for iid, i in self.c.inst.items():
            if "instance-id" in f and iid not in f["instance-id"]:
                continue
            if "client-token" in f and i["token"] not in f["client-token"]:
                continue
            if "tag:RunId" in f and (i["tags"]["RunId"] not in f["tag:RunId"]
                                     or i["state"] not in f["instance-state-name"]):
                continue
            d = {"InstanceId": iid, "State": {"Name": i["state"]}, "Placement": {"AvailabilityZone": i["az"]},
                 "Tags": [{"Key": k, "Value": v} for k, v in i["tags"].items()]}
            if i["reason"]:
                d["StateReason"] = {"Code": i["reason"]}
            out.append(d)
        return {"Reservations": [{"Instances": out}]}

    def terminate_instances(self, InstanceIds):
        for iid in InstanceIds:
            self.c.inst[iid]["state"] = "terminated"
        return {}


class _FakeS3:
    def __init__(self, c):
        self.c = c

    def get_object(self, Bucket, Key):
        self.c.tick()
        if Key not in self.c.objects:
            raise cerr("NoSuchKey", "GetObject")
        return {"Body": io.BytesIO(self.c.objects[Key])}

    def list_objects_v2(self, Bucket, Prefix, **kw):
        return {"IsTruncated": False, "Contents": [{"Key": k} for k in sorted(self.c.objects) if k.startswith(Prefix)]}


class _FakeSsm:
    def get_parameter(self, Name):
        return {"Parameter": {"Value": "ami-0fixture000000000"}}


def fake_runner(tmp_path, cloud, clk, **kw):
    return R.Runner(make_cfg(tmp_path, **kw), cloud.ec2, cloud.s3, cloud.ssm, clock=clk, sleep=clk.sleep,
                    out=lambda m: None)


def test_queue_runs_batch_under_vcpu_quota(tmp_path):
    """10 частей × 64 vCPU при квоте 384: одновременно не больше 6, остальные — по мере освобождения;
    бюджет проверяется перед каждым запуском; все части в самой дешёвой зоне."""
    clk = Clock()
    cloud = FakeCloud(clk)
    r = fake_runner(tmp_path, cloud, clk, max_parallel=10, vcpu_quota=384)
    res = r.run([make_job(name=f"part{i:02d}") for i in range(10)])
    assert all(v["status"] == "ok" for v in res.values()) and len(res) == 10
    assert cloud.max_live == 6 and cloud.max_vcpu == 384
    recs = ledger(tmp_path)
    launches = [x for x in recs if x["ev"] == "launch"]
    assert len(launches) == 10 and sum(1 for x in launches if x["t_launch"] > T0) == 4   # 4 — после освобождения
    for L in launches:                                          # перед каждым запуском — проверка бюджета этой части
        i = recs.index(L)
        assert recs[i - 1]["ev"] == "budget" and recs[i - 1]["job"] == L["job"] and recs[i - 1]["ok"]
    assert {x["az"] for x in launches} == {"us-east-1c"} and all(v["cheapest_zone"] for v in res.values())
    assert R.Ledger(tmp_path / "ledger.jsonl").summary()["open"] == []


def test_queue_waits_on_quota_error_without_own_machines(tmp_path):
    clk = Clock()
    cloud = FakeCloud(clk, fail=["VcpuLimitExceeded", "VcpuLimitExceeded"])   # квоту держат чужие машины
    r = fake_runner(tmp_path, cloud, clk)
    res = r.run([make_job(name="solo")])
    assert res["solo"]["status"] == "ok"
    assert [x["code"] for x in ledger(tmp_path, "launch_error")] == ["VcpuLimitExceeded"] * 2
    assert ledger(tmp_path, "launch")[0]["t_launch"] == T0 + 120             # две паузы по poll_s
    # предел ожидания исчерпан — задача снимается со статусом quota
    clk2 = Clock()
    cloud2 = FakeCloud(clk2, fail=["VcpuLimitExceeded"] * 5)
    r2 = fake_runner(tmp_path / "b", cloud2, clk2, quota_wait_polls=2)
    assert r2.run([make_job(name="solo")])["solo"]["status"] == "quota"


def test_queue_budget_wait_without_spamming_checks(tmp_path):
    clk = Clock()
    cloud = FakeCloud(clk)
    r = fake_runner(tmp_path, cloud, clk, ceiling=2.0)
    res = r.run([make_job(name="a"), make_job(name="b")])
    assert res["a"]["status"] == "ok" and res["b"]["status"] == "ok"
    budget = [(x["job"], x["itype"], x["ok"]) for x in ledger(tmp_path, "budget")]
    assert budget == [("a", "c7a.16xlarge", True), ("b", "c7a.16xlarge", False), ("b", "c7i.16xlarge", False),
                      ("b", "c7a.16xlarge", True)]                    # пересчёт — только после закрытия «a»
    la, lb = ledger(tmp_path, "launch")
    end_a = next(x for x in ledger(tmp_path, "end") if x["job"] == "a")
    assert lb["t_launch"] >= end_a["t_end"]


@pytest.mark.parametrize("reason", [R.SPOT_REASON, ""])        # "" — код причины не пришёл: решает метка 143
def test_spot_interruption_retry_without_duplicate_results(tmp_path, reason):
    """Прерывание: метка 143 первой попытки — не итог; повтор — в префикс _attempt2 и в другую зону;
    итог задачи — только из второй попытки, частичные файлы первой — отдельно в _lost/a1."""
    clk = Clock()
    cloud = FakeCloud(clk, interrupt={("a", 1): 120}, irq_reason=reason)
    r = fake_runner(tmp_path, cloud, clk)
    res = r.run([make_job(name="a")])
    a = res["a"]
    assert a["status"] == "ok" and a["attempts"] == 2 and "interrupted" not in a
    assert (a["az"], a["az_rank"], a["az_cheapest"], a["cheapest_zone"]) == ("us-east-1b", 1, "us-east-1c", False)
    run_dir = tmp_path / "runs" / RUN
    assert (run_dir / "a" / "result.txt").read_text() == "a a2"
    assert not (run_dir / "a" / "partial.txt").exists()
    assert (run_dir / "_lost" / "a1" / "a" / "partial.txt").read_text() == "partial"
    assert sorted(k for k in cloud.objects if "/_attempt2/" in k) == [
        f"{PREFIX}/gtu/{RUN}/_attempt2/a/_DONE.json", f"{PREFIX}/gtu/{RUN}/_attempt2/a/result.txt"]
    assert [(x["iid"], x["outcome"]) for x in ledger(tmp_path, "end")] == [("i-0001", "lost"), ("i-0002", "done")]
    assert len(cloud.runs) == 2 and cloud.runs[1]["SubnetId"] == SUBNETS["us-east-1b"]


def test_marker_of_other_attempt_is_ignored(tmp_path):
    clk = Clock()
    cloud = FakeCloud(clk, run_s=10 ** 9)
    r = fake_runner(tmp_path, cloud, clk)
    job = make_job(name="a")
    st, act = r.launch(job, attempt=2)
    r.active[act.iid] = act
    cloud.objects[r.key_prefix("a", attempt=2) + "/_DONE.json"] = json.dumps(
        {"attempt": 1, "phase": "done", "exit_code": 0}).encode()          # чужая (старая) метка
    r._poll(deque())
    assert "a" not in r.results and act.iid in r.active


# ----------------------------------------------------------------------------- параллельность и задачи
def test_parallel_limit_and_vcpu_quota(tmp_path):
    def fake_launch(self, job, attempt=1):
        iid = f"i-{job.name}"
        return "launched", R.Active(job, attempt, iid, job.instance_type, "us-east-1a", T0, T0 + 3600)

    jobs = [make_job(name=f"j{i}") for i in range(8)]
    for quota, limit, expect in ((384, 6, 6), (128, 6, 2), (384, 3, 3)):
        r = R.Runner(make_cfg(tmp_path, vcpu_quota=quota, max_parallel=limit), out=lambda m: None)
        r.launch = fake_launch.__get__(r)
        pending = deque((j, 1) for j in jobs)
        r._launch_some(pending)
        assert len(r.active) == expect and len(pending) == 8 - expect


def test_max_vcpu_cli_alias():
    p = R.build_parser()
    assert p.parse_args(["run", "j.json", "--max-vcpu", "128"]).vcpu_quota == 128
    assert p.parse_args(["run", "j.json", "--vcpu-quota", "64"]).vcpu_quota == 64
    assert p.parse_args(["plan", "j.json", "--assume-price", "1"]).vcpu_quota == 384


def test_load_jobs_validation(tmp_path):
    f = tmp_path / "jobs.json"
    f.write_text(json.dumps({"defaults": {"max_min": 45, "pip": ["numpy"]},
                             "jobs": [{"name": "a", "command": "true"},
                                      {"name": "b", "command": "true", "instance_type": "c7i.8xlarge",
                                       "fallback_types": []}]}))
    jobs = R.load_jobs(f)
    assert [(j.name, j.max_min, j.vcpu, j.types, j.image) for j in jobs] == [
        ("a", 45, 64, ["c7a.16xlarge", "c7i.16xlarge"], "al2023"), ("b", 45, 32, ["c7i.8xlarge"], "al2023")]
    for bad in ([{"name": "a", "command": "true", "typo": 1}],
                [{"name": "a", "command": "true"}, {"name": "a", "command": "true"}],
                [{"name": "a b", "command": "true"}],
                [{"name": "a", "command": "true", "max_min": 5}],
                [{"name": "a", "command": "true", "instance_type": "x9.huge"}],
                [{"name": "a", "command": "true", "pip": ["numpy; rm -rf /"]}],
                [{"name": "a", "command": "true", "image": "ubuntu2404", "dnf": ["gcc"]}],
                [{"name": "a", "command": "true", "image": "al2023", "apt": ["gmsh"]}]):
        f.write_text(json.dumps(bad))
        with pytest.raises(R.RunnerError):
            R.load_jobs(f)


@pytest.mark.parametrize("bad", [
    {"python": 3.12},                                            # число, а не строка (3.10 стало бы 3.1)
    {"python": "python3.12"},
    {"requirements": "../secret.txt"},
    {"requirements": "/etc/passwd"},
    {"container": "Ubuntu:24.04;rm"},
    {"container": "ubuntu:24.04", "pip": ["numpy"]},             # пакеты — внутри контейнера
    {"container": "ubuntu:24.04", "python": "3.12"},
    {"container": "ubuntu:24.04", "image": "ubuntu2404"},        # контейнер — на al2023
    {"container_setup": ["apt-get install -y x"]},               # без container
    {"container": "ubuntu:24.04", "container_setup": ["GTU_CMD_EOF"]},
    {"container": "ubuntu:24.04", "container_setup": ["  "]},
])
def test_new_fields_validation(bad):
    with pytest.raises(R.RunnerError):
        make_job(**bad)


def test_new_fields_accepted():
    j = make_job(container="opencfd/openfoam-default:2406", container_setup=["true"])
    assert (j.mode, j.image) == ("container", "al2023")
    assert make_job(container="public.ecr.aws/docker/library/ubuntu:24.04").mode == "container"
    assert make_job(python="3.12.7", requirements="routes/gtu/platform/docker/requirements.txt").mode == "python"
    assert make_job(python="3.12", image="ubuntu2404").image == "ubuntu2404"


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
               {"ebs_usd_gb_month": -1.0}, {"ipv4_usd_h": -1.0}, {"price_margin": 0.5}, {"ceiling": 0.0},
               {"quota_wait_polls": -1}, {"vcpu_quota": 0}):
        with pytest.raises(R.RunnerError):
            R.validate_config(make_cfg(tmp_path, **kw))
    R.validate_config(make_cfg(tmp_path))


def test_example_jobs_file_is_valid():
    jobs = R.load_jobs(HERE / "jobs.example.json")
    assert {j.name for j in jobs} == {"smoke", "explore-249"}
    assert {j.image for j in jobs} == {"al2023"}


# ----------------------------------------------------------------------------- совместимость старых задач
LEGACY_EXAMPLE = {   # jobs.example.json до v2 (коммит 55d1e8c): Ubuntu, системный Python
    "defaults": {"instance_type": "c7a.16xlarge", "fallback_types": ["c7i.16xlarge"], "image": "ubuntu2404",
                 "disk_gb": 50, "max_min": 30},
    "jobs": [{"name": "smoke", "max_min": 20,
              "command": "nproc > \"$GTU_OUT/nproc.txt\"\npython -c 'import sys, boto3; print(sys.version, "
                         "boto3.__version__)' > \"$GTU_OUT/python.txt\"\ngit rev-parse HEAD > \"$GTU_OUT/commit.txt\""},
             {"name": "explore-249", "command": "python routes/gtu/gtu_explore.py \"$GTU_OUT/gtu_explore_results.json\"",
              "pip": ["cantera", "scipy", "numpy"], "max_min": 40}]}
LEGACY_APT_NO_IMAGE = [{"name": "su2-old", "command": "true", "apt": ["build-essential", "openmpi-bin"],
                        "pip": ["numpy"], "max_min": 60}]                 # раньше образ по умолчанию был Ubuntu


def test_legacy_jobs_still_load_and_plan(tmp_path, capsys):
    for i, data in enumerate((LEGACY_EXAMPLE, LEGACY_APT_NO_IMAGE)):
        f = tmp_path / f"legacy{i}.json"
        f.write_text(json.dumps(data))
        jobs = R.load_jobs(f)
        assert {j.image for j in jobs} == {"ubuntu2404"} and {j.mode for j in jobs} == {"system"}
        for j in jobs:
            ud = R.build_userdata(j, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket=BUCKET,
                                  key_prefix="p/gtu/r/j")
            assert "apt-get -o DPkg::Lock::Timeout=600 -q update || exit 10" in ud
            assert 'python3 -m venv "$G/venv" || exit 12' in ud
        assert R.main(["plan", str(f), "--assume-price", "1.5", "--run-id", f"legacy{i}",
                       "--ledger", str(tmp_path / "l.jsonl"), "--out", str(tmp_path / "runs")]) == 0
    out = capsys.readouterr().out
    assert "canonical/ubuntu/server/24.04" in out and "/dev/sda1" in out


def _repo_job_files():
    pats = ("jobs*.json", "*aws_job*.json", "*runner_jobs*.json")
    return sorted({p for pat in pats for p in (REPO / "routes" / "gtu").rglob(pat) if "/runs/" not in str(p)})


def test_all_repo_job_files_plan(tmp_path):
    files = _repo_job_files()
    assert HERE / "jobs.example.json" in files and HERE / "jobs.smoke2.json" in files
    for f in files:
        assert R.main(["plan", str(f), "--assume-price", "1.5", "--run-id", "repo-plan",
                       "--ledger", str(tmp_path / "l.jsonl"), "--out", str(tmp_path / "runs")]) == 0, f
    assert not (tmp_path / "l.jsonl").exists()


# ----------------------------------------------------------------------------- второй пробный прогон (пункт Е)
def test_smoke2_jobs_file():
    jobs = {j.name: j for j in R.load_jobs(HERE / "jobs.smoke2.json")}
    assert set(jobs) == {"smoke-py", "smoke-ccx"}
    py, ccx = jobs["smoke-py"], jobs["smoke-ccx"]
    assert (py.mode, py.image, py.python, py.pip) == ("python", "al2023", "3.12", ["cantera", "numpy", "scipy"])
    assert (ccx.mode, ccx.image, ccx.container) == ("container", "al2023", "ubuntu:24.04")
    assert any("calculix-ccx" in c for c in ccx.container_setup) and ccx.container_setup[0].startswith("apt-get update")
    for j in (py, ccx):
        assert "routes/gtu/platform/aws/smoke2.py" in j.command and j.max_min <= 30
        ud = R.build_userdata(j, run_id=RUN, attempt=1, repo=R.DEFAULT_REPO, branch="main", bucket=BUCKET,
                              key_prefix=f"{PREFIX}/gtu/{RUN}/{j.name}", az="us-east-1c", az_cheapest="us-east-1c")
        p = subprocess.run(["bash", "-n"], input=ud, text=True)
        assert p.returncode == 0 and len(ud.encode()) < R.USERDATA_LIMIT
    for rel in ("routes/gtu/platform/aws/smoke2.py", "routes/gtu/gtu_cycle.py", "routes/gtu/platform/fea/disk.py"):
        assert (REPO / rel).exists(), rel


def test_plan_shows_zone_order_and_expected_price(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("ALPHA_O_S3", f"s3://{BUCKET}/{PREFIX}")
    rc = R.main(["plan", str(HERE / "jobs.smoke2.json"), "--assume-price", "0.955",
                 "--assume-zone-prices", "us-east-1a=0.955,us-east-1b=0.748,us-east-1c=0.706",
                 "--run-id", "plan-zones", "--ledger", str(tmp_path / "l.jsonl"), "--out", str(tmp_path / "runs")])
    out = capsys.readouterr().out
    assert rc == 0
    assert out.count("us-east-1c 0.7060 → us-east-1b 0.7480 → us-east-1a 0.9550") == 2
    assert "ожидается us-east-1c по 0.7060" in out and "MaxPrice 1.1460" in out
    assert '"SubnetId": "<подсеть по умолчанию us-east-1c>"' in out
    assert "Python 3.12 (uv)" in out and "контейнер ubuntu:24.04" in out
    new = 2 * (1.146 + 30 * 0.08 / 730 + 0.005) * (25 + 10) / 60          # две задачи по максимуму цены
    assert f"прогноз пакета {new:.2f}" in out
    assert BUCKET not in out


def test_smoke2_az_check_logic(tmp_path):
    def run(env, mode="ccx"):
        e = {"PATH": os.environ["PATH"], "GTU_OUT": str(tmp_path / mode), "NPROC": "2", **env}
        p = subprocess.run([sys.executable, str(HERE / "smoke2.py"), mode], env=e, capture_output=True, text=True,
                           timeout=300, cwd=tmp_path)
        return p.returncode, json.loads((tmp_path / mode / "az_check.json").read_text())

    if not shutil.which("ccx"):
        pytest.skip("ccx нет в контейнере")
    rc, az = run({"GTU_AZ": "us-east-1c", "GTU_AZ_CHEAPEST": "us-east-1c", "GTU_AZ_MACHINE": "us-east-1c"})
    assert rc == 0 and az["ok"]
    smoke = json.loads((tmp_path / "ccx" / "smoke_ccx.json").read_text())
    assert "Version" in smoke["ccx_version"] and smoke["disk"]["burst_margin"] > 1
    rc, az = run({"GTU_AZ": "us-east-1b", "GTU_AZ_CHEAPEST": "us-east-1c", "GTU_AZ_MACHINE": "us-east-1b"})
    assert rc == 3 and not az["ok"] and "не самую дешёвую" in az["why"]


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
    jobs = str(HERE / "jobs.example.json")
    rc = R.main(["plan", jobs, "--assume-price", "1.5", "--show-userdata", "--run-id", "cli-plan",
                 "--ledger", str(tmp_path / "ledger.jsonl"), "--out", str(tmp_path / "runs")])
    out = capsys.readouterr().out
    assert rc == 0 and "shutdown -P +20" in out and '"MarketType": "spot"' in out
    assert BUCKET not in out and FAKE_SECRET not in out and FAKE_AKID not in out
    assert "прогноз пакета 2.41" in out                       # (1.8 + диск) × (20+10 и 40+10 мин)
    assert "цены по зонам неизвестны" in out                  # без --assume-zone-prices зону выберет EC2
    assert not (tmp_path / "ledger.jsonl").exists()           # план не пишет в общий журнал
    write_ledger(tmp_path / "ledger.jsonl", HIST)
    assert R.main(["status", "--ledger", str(tmp_path / "ledger.jsonl")]) == 0
    out = capsys.readouterr().out
    assert "израсходовано ≈ 1.01" in out and "по ним 2.34" in out and "остаток 1996.64" in out   # 2.345 → 2.34 (float); потолок по умолчанию 2000 (решение 14.11.2026)


def test_cli_run_refuses_without_env_before_any_aws_call(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("ALPHA_O_S3", raising=False)
    monkeypatch.delenv("ALPHA_O_INSTANCE_PROFILE", raising=False)
    jobs = str(HERE / "jobs.example.json")
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
    add_subnets(aws)
    exp = dict(expected_run_params(), DryRun=True)
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
    assert f["ceiling"] == 2000.0 and f["report_step"] == 500.0 and f["vcpu_quota"] == 384
