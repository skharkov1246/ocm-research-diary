# -*- coding: utf-8 -*-
"""
NaCN 2500 t/y (Andrussow) modular plant, Tarkwa/Ahafo, Ghana.
Node: LOGISTICS MINIMISATION.
Package inventory: mass, envelope, China FOB cost, localisation eligibility.

ALL cost rates flagged SRC= are ASSUMPTIONS pending vendor quotation
(web-search budget exhausted this session; no primary source obtained).
Engineering quantities are derived from process mass balance (see balance.py).
"""

# ---- container geometry (ISO 668 / manufacturer data, standard values) ----
CONT = {
    # name: (int_len, int_wid, int_hgt m, int_vol m3, payload_t, tare_t)
    "40HC": (12.032, 2.352, 2.698, 76.3, 26.58, 3.90),
    "40GP": (12.032, 2.352, 2.393, 67.7, 26.68, 3.80),
    "20GP": ( 5.898, 2.352, 2.393, 33.2, 28.13, 2.25),
    "40OT": (12.030, 2.348, 2.330, 65.8, 26.63, 3.85),   # open top, roof removable
    "40FR": (12.060, 2.398, 2.140, 61.9, 39.50, 5.50),   # flat rack, collapsible
}
DIAG_40HC = (12.032**2 + 2.352**2 + 2.698**2) ** 0.5

# ---- freight rate assumptions (USD, all-in ocean incl. BAF/THC origin) ----
# SRC=ASSUMPTION: Far-East -> Tema CFR, mid-2026 band. NEEDS QUOTE.
RATE = {
    "ocean_40HC": 4000.0,      # band 3000-5500
    "ocean_40GP": 3800.0,
    "ocean_20GP": 2600.0,
    "ocean_40FR_ingauge": 6000.0,   # 1.5x
    "ocean_40FR_OOG": 12000.0,      # slot-loss billing, 3x  (band 2.5-4x)
    "bb_per_freight_ton": 220.0,    # breakbulk, 1t or 1m3 whichever greater
    "dest_port_40": 1100.0,    # Tema THC+GPHA+DO+agency+customs clearance/box
    "dest_port_OOG": 2600.0,   # + mobile crane, OOG handling, longer dwell
    "inland_40_tarkwa": 1200.0,  # Tema->Tarkwa ~330 km, 40' haulage
    "inland_OOG_tarkwa": 8500.0, # abnormal load: permit, survey, escort, utility lifts
    "insurance_pct": 0.011,    # marine cargo all-risk on CIF
}

# ---- Ghana local fabrication rates ----
# SRC=ASSUMPTION (Tema/Takoradi fabricator band). NEEDS QUOTE.
GH = {
    "struct_steel_fab_erect_per_t": 2700.0,  # band 2200-3200 $/t supply+fab+galv
    "cs_tank_shopfab_per_t":        3400.0,  # band 2900-4000 $/t CS tank, field-erected
    "lowclass_pipe_per_t":          4200.0,  # CW/service water/air, CS, sch40, incl fittings+install
    "cable_tray_per_t":             3600.0,
    "lv_cable_index":               1.10,    # local cable price vs China FOB (Tropical Cable/Nexans Tema)
    "concrete_per_m3":              175.0,   # ready-mix delivered, band 140-210
    "hdpe_liner_per_m2":            28.0,
}
# ---- China FOB reference rates for the same commodities ----
CN = {
    "struct_steel_per_t": 1600.0,
    "cs_tank_per_t":      2300.0,
    "lowclass_pipe_per_t":2500.0,
    "cable_tray_per_t":   2100.0,
}
