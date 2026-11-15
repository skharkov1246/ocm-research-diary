"""Воронка двигателя-заменителя, уровень 0: монтаж и ремонт «хоть уборщиком» -> балл 0-10.

Вход - перечень модулей замены (сменных блоков) концепции: масса, стыки (болтовой фланец или быстроразъёмный хомут),
трубки и кабели (с защитой от неправильного соединения или без), нужна ли центровка и балансировка, специнструмент,
частая ли это замена и какие модули надо снять, чтобы до него добраться. Перечень строит funnel_space.derive.

Балл (README узла, раздел «service_model»):
    S = 10 * (0,30*s_mass + 0,25*s_steps + 0,15*s_tools + 0,15*s_err + 0,15*s_diag)
    s_mass  = 1 - log10(M/25) / log10(2000/25), в пределах 0..1; M - масса самого тяжёлого модуля, кг
              (25 кг - подъём одним человеком без механизации; 2000 кг и тяжелее - 0);
    s_steps = 1 / (1 + шаги/10); шаги - среднее по частым заменам, вместе с модулями, снимаемыми для доступа;
    s_tools = 1 - 0,25 * число видов специнструмента (не меньше 0);
    s_err   = 1 / (1 + n/4); n - соединения частых замен, которые можно собрать неправильно;
    s_diag  = 0,5 + 0,25*[магнитные подшипники: встроенный контроль ротора] + 0,25*[модули с одним кодированным разъёмом].
Шаги модуля: болтовой стык - 1 + болты/10, хомут - 1; трубка - 1 (кодированный быстроразъём - 0,5); кабели - по одному
на разъём или 1 на общий кодированный разъём; центровка - 5; балансировка на месте - 5; строповка модуля тяжелее 25 кг - 2.
Все коэффициенты - оценка (README).

v3 (решение владельца 15.11.2026, п. 3): у модуля - где его ремонтируют (repair): «на месте» (модуль чинят или меняют на объекте),
«замена на месте, ремонт на заводе» (картридж газогенератора, силовая турбина), «на заводе» (оболочка: двигатель меняется целиком,
разборка - на заводе). Место ремонта в балл отдельным членом не входит: оно действует через модуль замены - у оболочки частая
замена - двигатель целиком (масса, стыки, строповка, центровка одновального), - и пишется в разбор (repair_sites, on_site_share).
"""
import math

M_ONE_PERSON, M_ZERO = 25.0, 2000.0
REPAIR_ON_SITE, REPAIR_SWAP, REPAIR_FACTORY = "на месте", "замена на месте, ремонт на заводе", "на заводе"
W = dict(mass=0.30, steps=0.25, tools=0.15, err=0.15, diag=0.15)
STEP = dict(joint=1.0, bolts_per_step=10.0, clamp=1.0, pipe=1.0, pipe_keyed=0.5, align=5.0, balance=5.0, sling=2.0)


def module_steps(m):
    s = 0.0
    for j in m.get("joints", []):
        s += STEP["clamp"] if j["type"] == "хомут" else STEP["joint"] + j.get("n_bolts", 0) / STEP["bolts_per_step"]
    s += m.get("pipes", 0) * (STEP["pipe_keyed"] if m.get("pipes_keyed") else STEP["pipe"])
    c = m.get("cables", 0)
    s += (1.0 if c > 0 else 0.0) if m.get("cables_keyed") else c
    s += STEP["align"] * bool(m.get("align")) + STEP["balance"] * bool(m.get("balance"))
    s += STEP["sling"] * (m.get("mass_kg", 0.0) > M_ONE_PERSON)
    return s


def module_wrongable(m):
    n = sum(1 for j in m.get("joints", []) if j["type"] != "хомут")
    n += 0 if m.get("pipes_keyed") else m.get("pipes", 0)
    n += 0 if m.get("cables_keyed") else m.get("cables", 0)
    n += 1 if m.get("align") else 0
    return n


def score(sv):
    mods = {m["name"]: m for m in sv["modules"]}
    if not mods:
        raise ValueError("нет модулей")
    M = max(m.get("mass_kg", 0.0) for m in mods.values())
    s_mass = 1.0 - math.log10(max(M, 1e-9) / M_ONE_PERSON) / math.log10(M_ZERO / M_ONE_PERSON)
    s_mass = min(1.0, max(0.0, s_mass))
    freq = [m for m in mods.values() if m.get("frequent")] or list(mods.values())
    rows = []
    for m in freq:
        acc = [mods[a] for a in m.get("access", []) if a in mods]
        st = module_steps(m) + sum(module_steps(a) for a in acc)
        wr = module_wrongable(m) + sum(module_wrongable(a) for a in acc)
        rows.append(dict(module=m["name"], steps=round(st, 2), wrongable=wr, access=[a["name"] for a in acc]))
    steps = sum(r["steps"] for r in rows) / len(rows)
    wrong = sum(r["wrongable"] for r in rows) / len(rows)
    tools = set()
    for m in mods.values():
        tools.update(m.get("tools", []))
        if m.get("align"):
            tools.add("центровка (лазерный прибор, проставки)")
        if m.get("balance"):
            tools.add("балансировка на месте")
    s_steps = 1.0 / (1.0 + steps / 10.0)
    s_tools = max(0.0, 1.0 - 0.25 * len(tools))
    s_err = 1.0 / (1.0 + wrong / 4.0)
    s_diag = 0.5 + 0.25 * bool(sv.get("amb")) + 0.25 * bool(sv.get("single_connector"))
    S = 10.0 * (W["mass"] * s_mass + W["steps"] * s_steps + W["tools"] * s_tools + W["err"] * s_err + W["diag"] * s_diag)
    heavy = max(mods.values(), key=lambda m: m.get("mass_kg", 0.0))
    sites = {m["name"]: m.get("repair", REPAIR_ON_SITE) for m in mods.values()}
    on_site = sum(1 for m in freq if m.get("repair", REPAIR_ON_SITE) == REPAIR_ON_SITE) / len(freq)
    return dict(score=round(S, 3), s_mass=round(s_mass, 3), s_steps=round(s_steps, 3), s_tools=round(s_tools, 3),
                s_err=round(s_err, 3), s_diag=round(s_diag, 3), heaviest=dict(module=heavy["name"], mass_kg=round(M, 1)),
                steps=round(steps, 2), wrongable=round(wrong, 2), tools=sorted(tools), frequent=rows,
                repair_sites=sites, on_site_share=round(on_site, 3),
                modules=[dict(name=m["name"], mass_kg=round(m.get("mass_kg", 0.0), 1), steps=round(module_steps(m), 2),
                              repair=m.get("repair", REPAIR_ON_SITE)) for m in mods.values()])
