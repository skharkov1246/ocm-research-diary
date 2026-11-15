#!/usr/bin/env python3
"""Точка входа образа платформы ГТУ: запуск задания по JSON (AWS Batch или локально).

    entrypoint.py --job <путь | s3://…/inputs/job.json> [--index N] [--local-root DIR] [--dry-run]

Задание (job.json, пишет render.py; схема gtu-batch-job/1):
    run_id        ключ прогона
    run_prefix    s3://<бакет>/<префикс>/gtu/<run_id> — корень прогона (inputs/, outputs/, logs/, manifest.json)
    plan          путь плана от корня прогона (inputs/plan.csv): строка = точка; столбцы index, row_id, line, …
    mode          point — задание с индексом N считает строку index == N;
                  line  — все строки line == N подряд (линия частоты: рестарт от предыдущей точки)
    command       bash-команда; {поле} — из строки плана; запускается из $GTU_HOME (корень кода в образе)
    np            число процессов MPI ($GTU_NP)
    time_limit_s  свой лимит на всё задание (тайм-аут Batch минус запас на выгрузку)

Индекс: --index, иначе AWS_BATCH_JOB_ARRAY_INDEX, иначе 0. Попытка: AWS_BATCH_JOB_ATTEMPT (иначе 1).
Для каждой строки: каталог /work/<row_id>/out ($GTU_OUT), журнал команды, _status.json; выгрузка сразу после
точки в outputs/<row_id>/attempt-<n>/ и logs/<row_id>/attempt-<n>.log (метка _status.json — последней).
Прерывание (SIGTERM: тайм-аут Batch или прерывание спота; через 30 с придёт SIGKILL) — команда гасится, что есть
выгружается со статусом terminated, выход 143. В режиме line линия останавливается на первой неудачной точке
(граница устойчивости стационарного расчёта или сбой): строки после неё — skipped.
Код выхода: 0 — все строки ok, или линия дала хотя бы одну точку (partial); иначе код первой неудачи (≥1);
2 — ошибка задания (нет плана, нет строк для индекса, неверный JSON).

--local-root DIR — хранилище вместо S3: s3://bucket/key ↔ DIR/bucket/key (локальный прогон и тесты, без boto3).
--dry-run — показать строки и команды, ничего не запускать и не выгружать.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import os
import pathlib
import shutil
import signal
import subprocess
import sys
import threading
import time

SCHEMA = "gtu-batch-job/1"
EXIT_JOB_ERROR = 2
EXIT_TERMINATED = 143


# ---------------------------------------------------------------- хранилище
def split_s3(uri: str) -> tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"не s3-адрес: {uri}")
    bucket, _, key = uri[5:].partition("/")
    return bucket, key


class LocalStorage:
    """s3://bucket/key ↔ <root>/bucket/key — для локального прогона и тестов."""

    def __init__(self, root: str):
        self.root = pathlib.Path(root)

    def _p(self, uri: str) -> pathlib.Path:
        b, k = split_s3(uri)
        return self.root / b / k

    def get_text(self, uri: str) -> str:
        return self._p(uri).read_text(encoding="utf-8")

    def upload_file(self, path: pathlib.Path, uri: str) -> None:
        dst = self._p(uri)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dst)

    def download_prefix(self, uri: str, dest: pathlib.Path) -> int:
        src = self._p(uri.rstrip("/"))
        n = 0
        if src.is_dir():
            for f in src.rglob("*"):
                if f.is_file():
                    t = dest / f.relative_to(src)
                    t.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(f, t)
                    n += 1
        return n


class S3Storage:
    """S3 через boto3 под ролью задания (jobRoleArn); ключей в образе и задании нет."""

    def __init__(self):
        import boto3  # noqa: PLC0415 — нужен только в облаке
        self.s3 = boto3.client("s3")

    def get_text(self, uri: str) -> str:
        b, k = split_s3(uri)
        return self.s3.get_object(Bucket=b, Key=k)["Body"].read().decode("utf-8")

    def upload_file(self, path: pathlib.Path, uri: str) -> None:
        b, k = split_s3(uri)
        self.s3.upload_file(str(path), b, k)

    def download_prefix(self, uri: str, dest: pathlib.Path) -> int:
        b, k = split_s3(uri.rstrip("/") + "/")
        n = 0
        for page in self.s3.get_paginator("list_objects_v2").paginate(Bucket=b, Prefix=k):
            for obj in page.get("Contents", []):
                rel = obj["Key"][len(k):]
                if not rel or rel.endswith("/"):
                    continue
                t = dest / rel
                t.parent.mkdir(parents=True, exist_ok=True)
                self.s3.download_file(b, obj["Key"], str(t))
                n += 1
        return n


# ---------------------------------------------------------------- задание и план
def load_job(ref: str, store) -> dict:
    text = store.get_text(ref) if ref.startswith("s3://") else pathlib.Path(ref).read_text(encoding="utf-8")
    job = json.loads(text)
    if job.get("schema") != SCHEMA:
        raise ValueError(f"схема задания {job.get('schema')!r}, ожидалась {SCHEMA}")
    for k in ("run_id", "run_prefix", "plan", "mode", "command", "np", "time_limit_s"):
        if k not in job:
            raise ValueError(f"в задании нет поля {k}")
    if job["mode"] not in ("point", "line"):
        raise ValueError("mode: point или line")
    if not str(job["run_prefix"]).startswith("s3://"):
        raise ValueError("run_prefix должен быть s3://…")
    return job


def read_plan(text: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(text)))
    for need in ("index", "row_id", "line"):
        if rows and need not in rows[0]:
            raise ValueError(f"в плане нет столбца {need}")
    return rows


def select_rows(rows: list[dict], mode: str, index: int) -> list[dict]:
    if mode == "point":
        sel = [r for r in rows if int(r["index"]) == index]
    else:
        sel = sorted((r for r in rows if int(r["line"]) == index), key=lambda r: int(r["index"]))
    return sel


def render_command(template: str, row: dict) -> str:
    """{поле} — из строки плана; ${…} bash не трогается (поля плана — только [A-Za-z0-9_])."""
    out, i = [], 0
    while i < len(template):
        c = template[i]
        if c == "{" and (i == 0 or template[i - 1] != "$"):
            j = template.index("}", i)
            name = template[i + 1:j]
            if name not in row:
                raise KeyError(f"в строке плана нет поля {name}")
            out.append(str(row[name]))
            i = j + 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def newest_restart(out_dir: pathlib.Path) -> str:
    files = sorted(out_dir.rglob("restart_*.dat"), key=lambda p: p.stat().st_mtime) if out_dir.is_dir() else []
    return str(files[-1]) if files else ""


# ---------------------------------------------------------------- запуск
class Runner:
    def __init__(self, job: dict, store, index: int, attempt: int, work: pathlib.Path, home: str,
                 clock=time.monotonic):
        self.job, self.store, self.index, self.attempt = job, store, index, attempt
        self.work, self.home, self.clock = work, home, clock
        self.t0 = clock()
        self.proc: subprocess.Popen | None = None
        self.terminated = False

    def on_term(self, signum, _frame):
        self.terminated = True
        if self.proc and self.proc.poll() is None:
            try:
                os.killpg(self.proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    @staticmethod
    def _expire(proc: subprocess.Popen, flag: threading.Event, grace_s: float = 20.0) -> None:
        """Срок вышел: SIGTERM группе процессов, через grace_s — SIGKILL (SU2 успевает дописать рестарт)."""
        if proc.poll() is not None:
            return
        flag.set()
        try:
            os.killpg(proc.pid, signal.SIGTERM)
            proc.wait(timeout=grace_s)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    def time_left(self) -> int:
        return int(self.job["time_limit_s"] - (self.clock() - self.t0))

    def uri(self, *parts: str) -> str:
        return "/".join([self.job["run_prefix"].rstrip("/"), *parts])

    def run_row(self, row: dict, prev: dict | None, prev_out: pathlib.Path | None) -> dict:
        rid = row["row_id"]
        rdir = self.work / rid
        out = rdir / "out"
        out.mkdir(parents=True, exist_ok=True)
        log = rdir / f"attempt-{self.attempt}.log"
        cmd = render_command(self.job["command"], row)
        env = dict(os.environ)
        left = self.time_left()
        env.update({
            "GTU_RUN_ID": self.job["run_id"], "GTU_ROW_ID": rid, "GTU_OUT": str(out), "GTU_NP": str(self.job["np"]),
            "GTU_ATTEMPT": str(self.attempt), "GTU_INPUTS": str(self.work / "inputs"),
            "GTU_TIME_LEFT": str(max(left - 60, 0)), "OMP_NUM_THREADS": "1",
            "GTU_PREV_OUT": str(prev_out) if prev_out else "",
            "GTU_PREV_RESTART": newest_restart(prev_out) if prev_out else "",
            "GTU_PREV_STATE": (prev or {}).get("state", "") if prev_out else "",
        })
        st = {"row_id": rid, "index": int(row["index"]), "line": int(row["line"]), "attempt": self.attempt,
              "run_id": self.job["run_id"], "command": cmd, "git_sha": os.environ.get("GTU_GIT_SHA", "unknown"),
              "batch_job_id": os.environ.get("AWS_BATCH_JOB_ID", ""), "batch_ce": os.environ.get("AWS_BATCH_CE_NAME", ""),
              "batch_jq": os.environ.get("AWS_BATCH_JQ_NAME", ""),
              "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        t1 = self.clock()
        if left <= 60:
            code, status = 124, "no_time"
        else:
            with open(log, "w", encoding="utf-8") as lf:
                lf.write(f"$ {cmd}\n")
                lf.flush()
                self.proc = subprocess.Popen(["bash", "-c", cmd], cwd=self.home, env=env, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT, start_new_session=True, text=True,
                                             errors="replace")
                timed_out = threading.Event()
                watchdog = threading.Timer(max(left - 30, 1), self._expire, args=(self.proc, timed_out))
                watchdog.daemon = True
                watchdog.start()
                for line in self.proc.stdout:          # в журнал и в stdout (CloudWatch Logs)
                    lf.write(line)
                    sys.stdout.write(line)
                code = self.proc.wait()
                watchdog.cancel()
            if self.terminated:
                status = "terminated"
            elif timed_out.is_set():
                status = "timeout"
            else:
                status = "ok" if code == 0 else "failed"
        st.update({"exit_code": code, "status": status, "seconds": round(self.clock() - t1, 1),
                   "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
        (out / "_status.json").write_text(json.dumps(st, ensure_ascii=False, indent=1), encoding="utf-8")
        self.upload(rid, out, log)
        return st

    def upload(self, rid: str, out: pathlib.Path, log: pathlib.Path) -> None:
        base = self.uri("outputs", rid, f"attempt-{self.attempt}")
        if log.exists():
            self.store.upload_file(log, self.uri("logs", rid, log.name))
        files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "_status.json")
        for p in files:
            self.store.upload_file(p, f"{base}/{p.relative_to(out).as_posix()}")
        self.store.upload_file(out / "_status.json", f"{base}/_status.json")   # метка — последней


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Запуск задания платформы ГТУ по JSON (AWS Batch / локально)")
    ap.add_argument("--job", required=True, help="job.json: путь или s3://…")
    ap.add_argument("--index", type=int, default=None, help="индекс (по умолчанию AWS_BATCH_JOB_ARRAY_INDEX или 0)")
    ap.add_argument("--local-root", default=None, help="папка вместо S3: s3://b/k ↔ DIR/b/k")
    ap.add_argument("--work", default=os.environ.get("GTU_WORK", "/work"), help="рабочий каталог")
    ap.add_argument("--dry-run", action="store_true", help="только показать строки и команды")
    a = ap.parse_args(argv)

    store = LocalStorage(a.local_root) if a.local_root else S3Storage()
    index = a.index if a.index is not None else int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX", "0"))
    attempt = int(os.environ.get("AWS_BATCH_JOB_ATTEMPT", "1"))
    try:
        job = load_job(a.job, store)
        rows = select_rows(read_plan(store.get_text(job["run_prefix"].rstrip("/") + "/" + job["plan"])),
                           job["mode"], index)
        if not rows:
            raise ValueError(f"в плане нет строк для индекса {index} (mode={job['mode']})")
        cmds = [render_command(job["command"], r) for r in rows]
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as e:
        print(f"[gtu] ошибка задания: {e}", file=sys.stderr)
        return EXIT_JOB_ERROR
    print(f"[gtu] прогон {job['run_id']} индекс {index} попытка {attempt} режим {job['mode']}: {len(rows)} строк")
    if a.dry_run:
        for r, c in zip(rows, cmds):
            print(f"[gtu] {r['row_id']}: {c}")
        return 0

    work = pathlib.Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    n_in = store.download_prefix(job["run_prefix"].rstrip("/") + "/inputs", work / "inputs")
    print(f"[gtu] входы: {n_in} файлов в {work / 'inputs'}")
    run = Runner(job, store, index, attempt, work, os.environ.get("GTU_HOME", os.getcwd()))
    signal.signal(signal.SIGTERM, run.on_term)

    results, prev, prev_out, first_bad = [], None, None, None
    for k, row in enumerate(rows):
        if run.terminated or (first_bad is not None and job["mode"] == "line"):
            results.append({"row_id": row["row_id"], "status": "skipped"})
            continue
        st = run.run_row(row, prev, prev_out)
        results.append(st)
        print(f"[gtu] {st['row_id']}: {st['status']} (код {st['exit_code']}, {st['seconds']} с)")
        if st["status"] == "ok":
            prev, prev_out = row, work / row["row_id"] / "out"
        elif first_bad is None:
            first_bad = st
    ok = sum(1 for r in results if r.get("status") == "ok")
    print(f"[gtu] итог: ok {ok} из {len(rows)}")
    if run.terminated:
        return EXIT_TERMINATED
    if first_bad is None or (job["mode"] == "line" and ok >= 1):
        return 0
    return int(first_bad["exit_code"]) or 1


if __name__ == "__main__":
    sys.exit(main())
