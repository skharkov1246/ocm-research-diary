# Узел bench: протокол сравнения CPU и GPU

Узел платформы трека «Простая литая ГТУ-заменитель» (issue #59), пункт 5 задания «Расчёт промышленной газовой
турбины на AWS»: по каким измерениям решать, переносить ли ядро SU2 на CUDA и когда переходить к нескольким
узлам. Сам протокол — `PROTOCOL.md`.

Состояние на 14.11.2026: прогонов по протоколу не было, в облаке ничего не запускалось. Готовы протокол,
конфигурация, разбор журналов SU2 с тестами и порог окупаемости по тарифам us-east-1.

| файл | что это |
|---|---|
| `PROTOCOL.md` | что умеют SU2, OpenFOAM, CalculiX на GPU; задача бенчмарка; что мерить; допуски точности; правило решения; переход к нескольким узлам |
| `bench.toml` | задача, критерии сходимости, допуски, правила решения, 8 конфигураций; тестовые значения помечены EXAMPLE, с владельцем данных и влиянием на точность (`[[example_data]]`) |
| `bench_parse.py` | журналы SU2 (`history*.csv`, `*.log`, `nvidia-smi`, `nsys`, `/usr/bin/time -v`) → `bench_summary.csv` и сверка с эталоном `bench_compare.csv` |
| `break_even.py` | `tariffs` (Price List AWS + Spot Instance Advisor), `mesh` (узлы и рёбра сеток), `compute` → `break_even.json` |
| `tariffs_us-east-1.json` | тарифы us-east-1 и eu-north-1 от 14.11.2026 с sha256 исходных файлов |
| `mesh_stats.json` | узлы и рёбра сеток узла geom (core130 coarse/medium, taurus60 coarse) |
| `break_even.json` | порог S* по парам GPU × CPU × рынок, доля ядра по Амдалу, оценка видеопамяти и обмена, упор в полосу памяти, сверка с цифрами задания |
| `test_bench.py` | 17 тестов на придуманных журналах (SYNTHETIC) |

```bash
python3 -m pytest -q test_bench.py
python3 break_even.py compute --bench bench.toml --out break_even.json
python3 bench_parse.py --bench bench.toml --runs runs --out bench_summary.csv --compare bench_compare.csv
```

Папки прогонов `runs/` в репозиторий не идут (`.gitignore`); сводки CSV и обновлённый `break_even.json` после
прогонов — идут.
