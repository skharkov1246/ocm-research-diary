"""Линейка типоразмеров двигателя-заменителя KV-1 … KV-n.

Решение владельца 14.11.2026: типоразмеры газогенератора называются KV-1, KV-2 и
так далее — порядковыми номерами по росту мощности на валу. Основа — мягкий допуск
пакета (решение того же дня) и охват 75 % расхода: frames["мягкий"] из
gtu_family_results.json. Номер типоразмера — его место по мощности на валу верхнего
(самого мощного) исполнения, которое он закрывает.

Скрипт ничего не пересчитывает в цикле: он читает результат gtu_family.py и
раскладывает исполнения по типоразмерам тем же правилом, что и отбор типоразмеров
(расход двигателя-заменителя в пределах [G_min; G_frame] с допуском на округление).

Запуск: python3 gtu_lineup.py → gtu_lineup_results.json
"""
import datetime
import json
import os
import re

os.chdir(os.path.dirname(os.path.abspath(__file__)))

SCENARIO = "мягкий"   # решение владельца 14.11.2026
TOL = 0.051           # расходы в результате округлены до 0,1 (типоразмер) и 0,01 (исполнение)


def base_model(m):
    return re.sub(r"\s*(\(.*?\))", "", m).strip()


def lineup(fam, scenario=SCENARIO):
    out = []
    for f in fam["frames"][scenario]:
        mem = [s for s in fam["slots"] if s.get(scenario) and f["G_min"] - TOL <= s[scenario]["G"] <= f["G_frame"] + TOL]
        if len(mem) != f["n"]:
            raise ValueError(f"типоразмер {f}: исполнений найдено {len(mem)}, в результате {f['n']}")
        out.append(dict(
            G_frame_kg_s=f["G_frame"], G_min_kg_s=f["G_min"], n=f["n"],
            P_shaft_MW_min=round(min(s["P_shaft_kW"] for s in mem) / 1000, 2),
            P_shaft_MW_max=round(max(s["P_shaft_kW"] for s in mem) / 1000, 2),
            members=[dict(model=s["model"], variant=s["variant"], P_kW=s["P_kW"], P_shaft_kW=s["P_shaft_kW"],
                          arch=s[scenario]["arch"], cool=s[scenario]["cool"], G_kg_s=s[scenario]["G"])
                     for s in sorted(mem, key=lambda s: s["P_shaft_kW"])],
            models=sorted({base_model(s["model"]) for s in mem}),
        ))
    out.sort(key=lambda x: (x["P_shaft_MW_max"], x["G_frame_kg_s"]))
    for i, x in enumerate(out, 1):
        x["name"] = f"KV-{i}"
    return out


if __name__ == "__main__":
    fam = json.load(open("gtu_family_results.json", encoding="utf-8"))
    res = dict(
        method="Раскладка исполнений по типоразмерам газогенератора из gtu_family_results.json (сценарий «мягкий», "
               "охват 75 % расхода); номер KV-n — место типоразмера по мощности на валу верхнего исполнения "
               "(решение владельца 14.11.2026).",
        software="python3",
        date=datetime.date.today().isoformat(),
        inputs=["gtu_family_results.json"],
        scenario=SCENARIO,
        lineup=lineup(fam),
    )
    json.dump(res, open("gtu_lineup_results.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    for x in res["lineup"]:
        print(f"{x['name']:5s} G {x['G_min_kg_s']:>6.1f}…{x['G_frame_kg_s']:<6.1f} кг/с  "
              f"P вала {x['P_shaft_MW_min']:>5.2f}…{x['P_shaft_MW_max']:<5.2f} МВт  n={x['n']:>2}  {', '.join(x['models'])}")
