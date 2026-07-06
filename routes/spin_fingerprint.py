#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/spin_fingerprint.py — скрин №1 off-label карты (METHANE_OFFLABEL_IDEAS.md):
«магнитный клапан» радикальных пар OCM. БЕЗ новых расчётов — майнинг спиновых
наблюдаемых из УЖЕ посчитанных JSON Этапов 15/16.

Физика. Магнитное поле управляет ветвлением радикальных реакций только там, где
конкурирующие каналы различаются спиновым характером у TS (радикально-парный
механизм, синглет/триплет-конверсия). Дешёвый отпечаток такой отзывчивости из
наших данных:
  1) Δ⟨S²⟩(TS) = spin-contamination расхождение каналов CH4 vs C2H6 вдоль
     HAT-профиля: рост ⟨S²⟩ над номиналом S(S+1) = появление радикальной пары
     (уходящий CH3•/C2H5• спин-скоррелирован с центром). Если каналы расходятся —
     поле бьёт по ним по-разному → клапан возможен.
  2) n_unpaired_beyond_spin и хвост NOON из CASSCF-profile — мультиреференсность
     у реперов (есть в emb-profile).
Выход: таблица + вердикт для статус-лога идей.
"""
import glob
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))


def load(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def geom_rows(path, sub):
    d = load(path)
    if not d:
        return []
    rows = []
    for i, p in enumerate(d.get("points", [])):
        rows.append({"sub": sub, "i": i, "rCH": p.get("rCH_pin"),
                     "s2": p.get("spin_square"), "e_h": p.get("e_h")})
    ts = d.get("ts_index")
    for r in rows:
        r["is_ts"] = (r["i"] == ts)
    return rows


def main():
    out = {"screen": "off-label #1 magnetic-valve spin fingerprint",
           "observable": "channel divergence of <S2> along HAT profile "
                         "(radical-pair character), CH4 vs C2H6",
           "sources": [], "tables": {}, "verdict": {}}
    S2_NOMINAL = 8.75  # септет 2S=5 -> S(S+1)=8.75 (Mn=O центр Этапов 15/16)

    for model, pref in (("stage16_embedded", "ocm_mnw_emb"),
                        ("stage15_bare", "ocm_mnw")):
        tab = {}
        for sub in ("ch4", "c2h6"):
            p = os.path.join(DIR, f"{pref}_{sub}_geom.json")
            rows = geom_rows(p, sub)
            if rows and any(r["s2"] for r in rows):
                tab[sub] = rows
                out["sources"].append(os.path.basename(p))
        if not tab:
            continue
        # контаминация = s2 - номинал; интересует max по профилю и на TS
        rep = {}
        for sub, rows in tab.items():
            cont = [(r["i"], round((r["s2"] or S2_NOMINAL) - S2_NOMINAL, 3),
                     bool(r.get("is_ts"))) for r in rows if r["s2"]]
            mx = max(cont, key=lambda t: t[1]) if cont else None
            ts = [c for c in cont if c[2]]
            rep[sub] = {"contam_per_point": cont, "max_contam": mx,
                        "ts_contam": ts[0] if ts else None}
        out["tables"][model] = rep
        if "ch4" in rep and "c2h6" in rep and rep["ch4"]["max_contam"] and \
           rep["c2h6"]["max_contam"]:
            dmax = round(rep["c2h6"]["max_contam"][1] - rep["ch4"]["max_contam"][1], 3)
            out["verdict"][model] = {
                "delta_max_contam_c2h6_minus_ch4": dmax,
                "reading": ("CHANNELS DIVERGE: C2H6 channel develops far more "
                            "radical-pair character -> magnetic lever plausible"
                            if dmax > 0.3 else
                            "channels similar -> weak magnetic leverage")}

    # NOON/мультиреференсность реперов из emb-profile (есть p0)
    for sub in ("ch4", "c2h6"):
        d = load(os.path.join(DIR, f"ocm_mnw_emb_{sub}_profile.json"))
        if d and d.get("points"):
            p0 = d["points"][0]
            out.setdefault("profile_multireference", {})[sub] = {
                "noon_tail": [n for n in p0.get("noon", []) if 0.05 < n < 1.95],
                "n_unpaired_beyond_spin": p0.get("n_unpaired_beyond_spin")}

    path = os.path.join(DIR, "spin_fingerprint_results.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
