# Генеративная геометрия тепла и потоков: «сосуды», решётки TPMS, оптимизация топологии

Узел трека «Простая литая ГТУ-заменитель» (`routes/gtu`, issue #59), папка `organic/thermal`. Задача владельца
14.11.2026: «футуристичная» геометрия, похожая на живой организм, — но как физика, а не украшение: сеть охлаждения
по законам ветвления, теплообменник на минимальных поверхностях, коллектор, выращенный оптимизацией топологии.
Приоритеты владельца: меньше деталей, Plug&Play в пакет Solar/Siemens, простой ремонт, затем КПД при малом размере.

Расчётная точка — KV-3 (Taurus 60 GS, `frame3/frame3_point_results.json`): π 13,5, T3 410,4 °C, T41 1000 °C,
25,55 кг/с, вал 14 950 об/мин. Нагрузка рекуператора — варианты «рекуператор ε + подогрев» из
`organic/cycles/cycles_results.json`. **Статус всего узла — оценка уровня облика.** Каждое число — из
results-JSON этой папки (поля `method`, `software`, `date` 2026-11-14, `inputs`, `assumptions`); внешние факты — с
источником и пометкой статуса; где что не получилось — сказано прямо (раздел «Что не получилось»).

## Главное (что показали расчёты)

ЗАПОЛНЯЕТСЯ ПО РЕЗУЛЬТАТАМ

## Файлы

| Файл | Что это |
|---|---|
| `voxflow.py` | воксельный решатель Стокса и теплообмена H1 в периодической ячейке (MAC, прямой SuperLU или MINRES + AMG) |
| `lbm.py` | решёточный Больцман D3Q19/D3Q7 (TRT): второй метод и течение с инерцией |
| `tpms_geom.py` | поля гироида и Schwarz-D, пористость, удельная поверхность (marching cubes), STL, упрощение сетки |
| `tpms_transport.py` | таблица f·Re и Nu_H1 по уровню поверхности (MAC; LBM — с инерцией), пул процессов |
| `tpms_recup.py` | рекуператор: противоток по участкам, переменная ячейка a(x) и смещение c(x), пластины для сравнения |
| `tpms.py` | сводка: геометрия → переносы → сверка методов → рекуператоры → STL, PNG → `tpms_results.json` |
| `vascular.py` | сети охлаждения: гидравлика (Ньютон) + 2D-стенка, лопатка СА и жаровая труба → `vascular_results.json` |
| `topo.py` | оптимизация топологии коллектора (Стокс — Бринкман, MMA, проекция) → `topo_results.json` |
| `make_jobs_thermal.py` → `jobs_thermal.json` | задачи AWS (формат `platform/aws/runner.py`) |
| `feed_thermal.py` → `feed_thermal.json` | блок `thermal` в `../feed_to_funnels.json` (оси и модели для воронок) |
| `test_thermal.py` | 10 проверок pytest на аналитике и балансах |
| `*.png`, `*.stl` | картинки и уменьшенная геометрия (полные STL — в scratchpad, см. «Большие файлы») |

## Как запускать

```bash
python3 -m venv --system-site-packages venv && venv/bin/pip install scikit-image trimesh pyamg fast-simplification
venv/bin/python -m pytest -q test_thermal.py                         # 10 passed, ~30 с
venv/bin/python tpms_transport.py --method mac --n 32 --shifts 0 0.5 --pe 0.01 10 20 --workers 3 --out tpms_transport_mac32.json
venv/bin/python tpms.py --mac tpms_transport_mac32.json --out . --big <scratchpad>
venv/bin/python vascular.py --out .
venv/bin/python topo.py --nx 80 --vol 0.30 --out .
python3 make_jobs_thermal.py && python3 ../../platform/aws/runner.py plan jobs_thermal.json --assume-price 0.955
python3 feed_thermal.py
```

Локально машина была общей (нагрузка 20–33 на 4 ядрах из-за соседних задач), поэтому время счёта в JSON —
при этой нагрузке. AWS не вызывался.

ЧАСТИ 1–3 — ЗАПОЛНЯЮТСЯ ПО РЕЗУЛЬТАТАМ
