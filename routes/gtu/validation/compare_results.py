#!/usr/bin/env python3
"""Сверка повторного расчёта с закоммиченным результатом (повторяемость эталонных расчётов трека ГТУ).

Два JSON одного скрипта (например, validation_fleet_results.json из ветки и тот же файл, пересчитанный на машине
AWS) обходятся параллельно: словари — по объединению ключей, списки — по индексу. Числа сравниваются с допуском
|новое − эталон| ≤ atol + rtol·|эталон|; строки и прочие значения — на равенство (расхождение строк считается, но
проверку не валит: в результатах бывают пути и версии). Служебные ключи (дата, версии ПО, время счёта, число
процессов, рабочие папки) пропускаются — список DEFAULT_IGNORE, дополняется --ignore.

Выход — небольшой JSON: число сравненных чисел, наибольшие абсолютное и относительное отклонения с путём,
первые расхождения, недостающие и лишние ключи, вердикт. Код выхода: 0 — в допуске, 3 — есть числа вне допуска
или разная структура, 2 — файл не прочитан.

  python3 compare_results.py NEW.json REF.json --out repro.json [--rtol 1e-6] [--atol 1e-9] [--ignore k1,k2]

Сбор итогов задач AWS (jobs_aws_validation.json, блок merge): --collect <runs/run_id> — все repro_*.json, timing.txt,
env.txt, commit.txt и последние строки pytest*.txt задач в одну сводку для ветки (без больших файлов).

  python3 compare_results.py --collect routes/gtu/platform/aws/runs/<run_id> --run-id <run_id> --out aws_repro_results.json
"""
from __future__ import annotations

import argparse
import datetime
import glob
import json
import math
import os
import sys

DEFAULT_IGNORE = {"date", "software", "runtime_s", "workers", "seconds", "wall_s", "workdir", "dir",
                  "created_utc", "generated_utc", "platform", "python", "host", "hostname"}


def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def compare(new, ref, rtol: float = 1e-6, atol: float = 1e-9, ignore=frozenset(DEFAULT_IGNORE), max_list: int = 20):
    st = dict(numbers=0, numbers_bad=0, strings=0, strings_diff=0, nonfinite_mismatch=0,
              missing_in_new=[], extra_in_new=[], type_mismatch=[], length_mismatch=[], bad=[], string_diffs=[],
              max_abs=dict(value=0.0, path=None), max_rel=dict(value=0.0, path=None))

    def note(lst, item):
        if len(lst) < max_list:
            lst.append(item)

    def walk(a, b, path):
        if isinstance(b, dict):
            if not isinstance(a, dict):
                note(st["type_mismatch"], path or "/")
                return
            for k in sorted(set(a) | set(b), key=str):
                if k in ignore:
                    continue
                p = f"{path}/{k}"
                if k not in a:
                    note(st["missing_in_new"], p)
                elif k not in b:
                    note(st["extra_in_new"], p)
                else:
                    walk(a[k], b[k], p)
            return
        if isinstance(b, list):
            if not isinstance(a, list):
                note(st["type_mismatch"], path)
                return
            if len(a) != len(b):
                note(st["length_mismatch"], dict(path=path, new=len(a), ref=len(b)))
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{path}[{i}]")
            return
        if _is_num(b) or b is None and _is_num(a):
            if not (_is_num(a) and _is_num(b)):
                note(st["type_mismatch"], path)
                return
            st["numbers"] += 1
            fa, fb = float(a), float(b)
            if not (math.isfinite(fa) and math.isfinite(fb)):
                if not (fa == fb or (math.isnan(fa) and math.isnan(fb))):
                    st["nonfinite_mismatch"] += 1
                    st["numbers_bad"] += 1
                    note(st["bad"], dict(path=path, new=a, ref=b))
                return
            d = abs(fa - fb)
            rel = d / abs(fb) if fb != 0 else (0.0 if d == 0 else math.inf)
            if d > st["max_abs"]["value"]:
                st["max_abs"] = dict(value=d, path=path, new=a, ref=b)
            if rel > st["max_rel"]["value"] and d > atol:
                st["max_rel"] = dict(value=rel, path=path, new=a, ref=b)
            if d > atol + rtol * abs(fb):
                st["numbers_bad"] += 1
                note(st["bad"], dict(path=path, new=a, ref=b, abs=d, rel=rel))
            return
        if isinstance(b, str) or isinstance(a, str):
            st["strings"] += 1
            if a != b:
                st["strings_diff"] += 1
                note(st["string_diffs"], dict(path=path, new=str(a)[:200], ref=str(b)[:200]))
            return
        if a != b:  # bool, None
            note(st["type_mismatch"], path)

    walk(new, ref, "")
    structure_ok = not (st["missing_in_new"] or st["extra_in_new"] or st["type_mismatch"] or st["length_mismatch"])
    st["structure_ok"] = structure_ok
    st["ok"] = structure_ok and st["numbers_bad"] == 0 and st["numbers"] > 0
    for k in ("max_abs", "max_rel"):
        v = st[k]["value"]
        st[k]["value"] = None if v == math.inf else v
    return st


def _tail(path: str, n: int = 3) -> list[str]:
    try:
        return [ln.rstrip() for ln in open(path, encoding="utf-8", errors="replace").read().splitlines() if ln.strip()][-n:]
    except OSError:
        return []


def collect(run_dir: str, run_id: str = "") -> dict:
    """Сводка задач прогона: сверки repro_*.json (без списков подробностей сверх первых), время, окружение, тесты."""
    jobs = {}
    for d in sorted(glob.glob(os.path.join(run_dir, "*"))):
        if not os.path.isdir(d) or os.path.basename(d).startswith("_"):
            continue
        reps = sorted(glob.glob(os.path.join(d, "repro_*.json")))
        if not reps:
            continue
        j = dict(repro={}, timing=_tail(os.path.join(d, "timing.txt"), 20), env=_tail(os.path.join(d, "env.txt"), 10),
                 commit=(_tail(os.path.join(d, "commit.txt"), 1) or [None])[0],
                 pytest={os.path.basename(p): _tail(p, 2) for p in sorted(glob.glob(os.path.join(d, "pytest*.txt")))})
        for r in reps:
            try:
                x = json.load(open(r, encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as e:
                j["repro"][os.path.basename(r)] = dict(ok=False, error=type(e).__name__)
                continue
            j["repro"][os.path.basename(r)] = {k: x.get(k) for k in (
                "label", "ok", "structure_ok", "numbers", "numbers_bad", "strings_diff", "max_abs", "max_rel", "bad",
                "missing_in_new", "extra_in_new", "length_mismatch", "rtol", "atol")}
        jobs[os.path.basename(d)] = j
    ok = bool(jobs) and all(r.get("ok") for j in jobs.values() for r in j["repro"].values())
    return dict(title="Повтор эталонных расчётов трека ГТУ на машинах AWS: сверка с веткой", run_id=run_id,
                method="compare_results.py на машине (повтор против закоммиченного JSON) + --collect по папке прогона",
                created_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                ok=ok, jobs=jobs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("new", nargs="?", help="пересчитанный JSON")
    ap.add_argument("ref", nargs="?", help="эталон (закоммиченный JSON)")
    ap.add_argument("--collect", default="", help="папка прогона runs/<run_id>: собрать repro_*.json задач")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--out", help="куда записать итог сверки (JSON); без него — печать")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--ignore", default="", help="доп. ключи через запятую (пропускаются на любом уровне)")
    ap.add_argument("--label", default="", help="подпись сверки")
    a = ap.parse_args(argv)
    if a.collect:
        doc = collect(a.collect, a.run_id)
        if not doc["jobs"]:
            print(f"[сверка] в {a.collect} нет repro_*.json", file=sys.stderr)
            return 2
        txt = json.dumps(doc, ensure_ascii=False, indent=1, default=str)
        if a.out:
            with open(a.out, "w", encoding="utf-8") as f:
                f.write(txt + "\n")
        for job, j in doc["jobs"].items():
            for name, r in j["repro"].items():
                print(f"[сверка] {job}/{name}: {'OK' if r.get('ok') else 'РАСХОЖДЕНИЕ'}, чисел {r.get('numbers')}, "
                      f"вне допуска {r.get('numbers_bad')}")
        return 0 if doc["ok"] else 3
    if not (a.new and a.ref):
        ap.error("нужны NEW и REF или --collect")
    try:
        new = json.load(open(a.new, encoding="utf-8"))
        ref = json.load(open(a.ref, encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        print(f"[сверка] файл не прочитан: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    ignore = set(DEFAULT_IGNORE) | {k.strip() for k in a.ignore.split(",") if k.strip()}
    st = compare(new, ref, a.rtol, a.atol, frozenset(ignore))
    doc = dict(label=a.label, new=a.new, ref=a.ref, rtol=a.rtol, atol=a.atol, ignored_keys=sorted(ignore),
               method="compare_results.py: параллельный обход двух JSON; числа |new-ref| <= atol + rtol*|ref|",
               created_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), **st)
    txt = json.dumps(doc, ensure_ascii=False, indent=1, default=str)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(txt + "\n")
    print(f"[сверка] {a.label or a.new}: чисел {st['numbers']}, вне допуска {st['numbers_bad']}, "
          f"макс. отн. {st['max_rel']['value']} ({st['max_rel']['path']}), строк различно {st['strings_diff']}, "
          f"структура {'ok' if st['structure_ok'] else 'РАЗНАЯ'} -> {'OK' if st['ok'] else 'РАСХОЖДЕНИЕ'}")
    return 0 if st["ok"] else 3


if __name__ == "__main__":
    sys.exit(main())
