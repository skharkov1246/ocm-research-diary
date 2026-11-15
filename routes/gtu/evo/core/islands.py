"""Детерминированные острова: сид острова, деление островов на части по номеру острова, прогон части, слияние архивов
с контролем полноты.

Правила:
  * сид острова i — sha256(base_seed | i | salt), первые 8 байт: не зависит от числа частей, процессов и машин;
  * остров i относится к части i mod n_parts — разбиение ровно одно, каждая часть знает свои острова без обмена;
  * острова внутри прогона не обмениваются особями (иначе итог зависел бы от деления на части); обмен — по эпохам:
    слитый архив эпохи e — начальные родители эпохи e + 1 (`seed_genomes`), так же входит выбор владельца;
  * слияние проверяет, что пришли все ожидаемые острова с одной и той же конфигурацией (хэш), и пишет отчёт:
    каких нет, какие лишние или повторные, сколько оценок; неполное слияние — ошибка, если не разрешено явно.
"""
import hashlib
import json
import multiprocessing as mp

import numpy as np

from . import archive as A


def island_seed(base_seed, island, salt=""):
    h = hashlib.sha256(f"{int(base_seed)}|{int(island)}|{salt}".encode("utf-8")).digest()
    return int.from_bytes(h[:8], "little")


def rng_for(base_seed, island, salt=""):
    return np.random.Generator(np.random.PCG64(island_seed(base_seed, island, salt)))


def islands_of_part(n_islands, n_parts, part):
    if not (1 <= n_parts <= n_islands):
        raise ValueError(f"частей {n_parts} при островах {n_islands}: нужно 1 ≤ частей ≤ островов")
    if not (0 <= part < n_parts):
        raise ValueError(f"часть {part} вне [0, {n_parts})")
    return [i for i in range(n_islands) if i % n_parts == part]


def config_hash(cfg):
    s = json.dumps(cfg, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def _call(args):
    fn, island, kw = args
    return fn(island, **kw)


def run_part(run_island, n_islands, n_parts, part, workers=1, **kw):
    """run_island(island, **kw) → dict (архив острова: archive.to_dict() + служебные поля). Порядок результатов —
    по номеру острова; от числа процессов не зависит."""
    ids = islands_of_part(n_islands, n_parts, part)
    jobs = [(run_island, i, kw) for i in ids]
    if workers <= 1 or len(ids) <= 1:
        return [_call(j) for j in jobs]
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
    with ctx.Pool(min(workers, len(ids))) as pool:
        return pool.map(_call, jobs, chunksize=1)


def merge_archives(island_dicts, expected_islands, cfg_hash=None, allow_partial=False):
    """Слияние архивов островов. island_dicts — список {island, config_hash, evals, archive: Archive.to_dict()}.
    Возврат (Archive, report). Правило ячейки — как в Archive.try_add (больший фитнес, при равенстве меньший uid),
    поэтому итог не зависит от порядка входов."""
    exp = sorted(int(i) for i in expected_islands)
    seen, dup, foreign, bad_hash = {}, [], [], []
    for d in island_dicts:
        i = int(d["island"])
        if cfg_hash is not None and d.get("config_hash") != cfg_hash:
            bad_hash.append(i)
            continue
        if i not in exp:
            foreign.append(i)
            continue
        if i in seen:
            dup.append(i)
            if _digest(seen[i]) != _digest(d):
                raise ValueError(f"остров {i} пришёл дважды с разным содержимым — разные прогоны в одной папке")
            continue
        seen[i] = d
    missing = [i for i in exp if i not in seen]
    report = dict(expected=len(exp), present=len(seen), missing=missing, duplicates=sorted(set(dup)), foreign=sorted(set(foreign)),
                  config_mismatch=sorted(set(bad_hash)), complete=not missing and not bad_hash,
                  evals=int(sum(int(d.get("evals", 0)) for d in seen.values())))
    if not report["complete"] and not allow_partial:
        raise ValueError(f"слияние неполное: нет островов {missing[:20]}{'…' if len(missing) > 20 else ''}, "
                         f"чужая конфигурация {bad_hash[:20]}")
    if not seen:
        raise ValueError("нечего сливать")
    grid = A.GridSpec.from_dict(next(iter(seen.values()))["archive"]["grid"])
    merged = A.Archive(grid)
    per = {}
    for i in sorted(seen):
        a = A.Archive.from_dict(seen[i]["archive"])
        if a.grid.to_dict() != grid.to_dict():
            raise ValueError(f"остров {i}: другая сетка дескрипторов")
        per[i] = dict(filled=len(a), best=(a.best()["fitness"] if len(a) else None))
        for c in a.cells():
            e = a.elites[c]
            merged.try_add(e["genome"], e["fitness"], e["desc"], e["objs"], e["meta"], e["uid"])
    report["per_island"] = per
    report["merged_filled"] = len(merged)
    report["merged_coverage"] = merged.coverage()
    report["status"] = "полное" if report["complete"] else f"НЕПОЛНОЕ: нет {len(missing)} островов из {len(exp)}"
    return merged, report


def _digest(d):
    return hashlib.sha256(json.dumps(d.get("archive", {}), sort_keys=True, default=str).encode()).hexdigest()
