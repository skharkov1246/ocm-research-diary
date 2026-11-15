#!/usr/bin/env python3
"""Рендер элит семейства «лопатки и колёса» и лист-галерея для выбора владельцем (интерактивная эволюция).

Идёт в окружении с pyvista (OSMesa, вне экрана); при её отсутствии — matplotlib. Модели расчёта не импортирует:
читает сетки элит (NPZ из evo_blades.py) и их записи (mini_results.json или elites_*.json).

Три вида на элиту, монохром, как эскизы воронки: сбоку (ось вала горизонтально), с торца (вдоль оси), аксонометрия
3/4 (азимут 55°, возвышение 30°). Подпись листа: «оценка уровня облика, не конструкторская документация».

  python render_blades.py --results DIR/mini_results.json --meshes DIR/meshes --out DIR/render [--thumb-dir REPO_DIR]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import core_local as C   # noqa: E402

TITLES = {"comp1": "колесо ступени 1 (двусторонний вход)", "comp2": "колесо ступени 2", "turb": "ротор турбины, 3 блиска"}
TECH_RU = {"tubercles": "бугорки", "splitter_offset": "смещённая разделительная", "serrations": "серрации", "tandem": "тандем",
           "riblets": "риблеты", "turb_bow": "составной навал", "turb_vortex": "закон закрутки"}


def load_meshes(p):
    z = np.load(p)
    n = int(z["n"][0])
    return [(z[f"V{i}"].astype(float), z[f"T{i}"].astype(int), float(z[f"g{i}"][0])) for i in range(n)]


def caption(comp, e):
    m = e["meta"]
    g = e["genes"]
    lines = []
    if comp == "turb":
        lines.append(f"η_p {m['eta']:.4f} (ур. 1 {m['eta_level1']:.4f}); m {m['mass_kg']:.0f} кг; J {m['J_kgm2']:.2f} кг·м²")
        lines.append(f"φ {g['phi']:.2f}, U {g['U']:.0f} м/с, zt {g['zt']}, сужение {g['taper']:.2f}, "
                     f"полость {m.get('hollow', 0):.0%}, n {g['vortex']:+.2f}")
        lines.append(f"навал {g['lean']:+.0f}°, составной {np.degrees(np.arctan(4 * g['bow'])):.0f}°, стрела {g['sweep']:+.0f}°, "
                     f"бугорки {g['tub_A']:.3f}; {'печать' if g['route'] == 0 else 'литьё'}")
    else:
        lines.append(f"η_p {m['eta']:.4f} (ур. 1 {m['eta_level1']:.4f}); запас помпажа {m['SM']:.0%}; m {m['mass_kg']:.0f} кг; "
                     f"J {m['J_kgm2']:.2f} кг·м²")
        sp = f"разд. {g['lsp']:.2f} @ {g['sp_off']:.2f} шага" if g["lsp"] > 0 else "без разделительных"
        lines.append(f"Z {g['Z']}, {sp}, β2b {m['beta2b']:.0f}°, навал {g['lean']:+.0f}°, стрела {g['sweep']:+.2f}")
        extra = []
        if g["tub_A"] > 0:
            extra.append(f"бугорки {g['tub_A']:.3f}×{g['tub_n']}")
        if g["ser_A"] > 0:
            extra.append(f"серрации {g['ser_A']:.2f}×{g['ser_n']}")
        if g["tand"] > 0:
            extra.append(f"тандем {g['tand_m']:.2f}/{g['tand_clk']:.2f}")
        if g.get("fish", 0) > 0:
            extra.append(f"«рыбья кость» {g['fish']:.2f}")
        lines.append((", ".join(extra) or "без приёмов") + f"; загрузка диска {m['u_disk']:.2f}; {'печать' if g['route'] == 0 else 'литьё'}")
    ef = e.get("effect", {})
    lines.append(f"новизна: в архиве {ef.get('novelty_archive', 0):.2f}, от типовых OEM {ef.get('novelty_oem', 0):.2f}")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--meshes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--thumb-dir", default=None)
    ap.add_argument("--size", type=int, default=420)
    ap.add_argument("--backend", default="auto")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    R = json.load(open(a.results, encoding="utf-8"))
    comps = R["comps"] if "comps" in R else {R["comp"]: dict(elites=R["elites"])}
    report = {}
    for comp, cd in comps.items():
        tiles = []
        eng = None
        for e in cd["elites"]:
            mp = os.path.join(a.meshes, f"{comp}_cell{e['cell']:04d}.npz")
            if not os.path.exists(mp):
                continue
            meshes = load_meshes(mp)
            paths, eng = C.render_views(meshes, os.path.join(a.out, f"{e['id']}_v{{}}.png"), (a.size, a.size), a.backend)
            lab = f"{e['id']}  (пригодность {e['f']:.2f})" if e["feasible"] else f"{e['id']}  НЕДОПУСТИМО: {e['meta'].get('why')}"
            tiles.append((lab, paths, caption(comp, e)))
        if not tiles:
            continue
        title = f"Семейство «лопатки и колёса» — {TITLES[comp]}: элиты мини-прогона (MAP-Elites → NSGA-II → CMA-ES)"
        footer = ("Оценка уровня облика, не конструкторская документация. Виды: сбоку (ось вала горизонтально), с торца, аксонометрия 3/4. "
                  "КПД — модель уровня 1 воронки + поправки приёмов (источники и статусы — README); прочность — быстрая оценка, "
                  "калиброванная по CalculiX. Выбор владельца: id в picks.json → evo_blades.py mini|islands --parents picks.json.")
        full = os.path.join(a.out, f"gallery_{comp}.png")
        C.gallery_sheet(tiles, full, title, footer, ncol=2, tile_px=(1050, 350))
        report[comp] = dict(sheet=full, n=len(tiles), engine=eng)
        if a.thumb_dir:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib import image as mpimg
            import matplotlib.pyplot as plt
            img = mpimg.imread(full)
            h, w = img.shape[:2]
            k = 1400 / w
            fig = plt.figure(figsize=(1400 / 100, h * k / 100), dpi=100)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.imshow(img)
            ax.set_axis_off()
            tp = os.path.join(a.thumb_dir, f"gallery_{comp}.png")          # уменьшенная копия (≤ 200 КБ после ужатия ниже)
            fig.savefig(tp, dpi=100)
            plt.close(fig)
            try:
                from PIL import Image
                im = Image.open(tp).convert("L")
                w_, h_ = im.size
                im = im.resize((1150, int(h_ * 1150 / w_)), Image.LANCZOS).quantize(colors=32).convert("L")
                im.save(tp, optimize=True)
            except ImportError:
                pass
            report[comp]["thumb"] = tp
            report[comp]["thumb_kb"] = round(os.path.getsize(tp) / 1024, 1)
    with open(os.path.join(a.out, "render_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    main()
