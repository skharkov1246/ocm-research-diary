"""Источники методов и приёмов модуля evo/core. Статус «по памяти» — ссылка библиографическая, первоисточник в
сессии не открывался; число из такого источника в расчёт не идёт без пометки «оценка». Свои числа трека — с файлом."""

SOURCES = {
    "map_elites": dict(ref="Mouret J.-B., Clune J. Illuminating search spaces by mapping elites. arXiv:1504.04909 (2015)",
                       status="по памяти"),
    "iso_line": dict(ref="Vassiliades V., Mouret J.-B. Discovering the elite hypervolume by leveraging interspecies correlation. GECCO 2018",
                     status="по памяти"),
    "nsga2": dict(ref="Deb K., Pratap A., Agarwal S., Meyarivan T. A fast and elitist multiobjective genetic algorithm: NSGA-II. "
                      "IEEE Trans. Evol. Comput. 6(2):182–197 (2002)", status="по памяти"),
    "zdt": dict(ref="Zitzler E., Deb K., Thiele L. Comparison of multiobjective evolutionary algorithms: empirical results. "
                    "Evol. Comput. 8(2):173–195 (2000)", status="по памяти; фронт ZDT1 f2 = 1 − √f1 проверен тестом"),
    "sbx_pm": dict(ref="Deb K., Agrawal R.B. Simulated binary crossover for continuous search space. Complex Systems 9 (1995); "
                       "Deb K. Multi-Objective Optimization using Evolutionary Algorithms. Wiley (2001)", status="по памяти"),
    "cmaes": dict(ref="Hansen N. The CMA Evolution Strategy: A Tutorial. arXiv:1604.00772 (2016)", status="по памяти; проверено тестом "
                      "(сфера, Розенброк)"),
    "novelty": dict(ref="Lehman J., Stanley K.O. Abandoning objectives: evolution through the search for novelty alone. "
                        "Evol. Comput. 19(2):189–223 (2011)", status="по памяти"),
    "ffd": dict(ref="Sederberg T.W., Parry S.R. Free-form deformation of solid geometric models. SIGGRAPH 1986", status="по памяти"),
    "smooth_min": dict(ref="полиномиальный smooth-min для SDF (И. Килес, заметки iquilezles.org «smooth minimum»)", status="по памяти"),
    "surface_nets": dict(ref="Gibson S.F.F. Constrained elastic surface nets. MICCAI 1998", status="по памяти"),
    "tpms": dict(ref="Schoen A.H. Infinite periodic minimal surfaces without self-intersections. NASA TN D-5541 (1970); "
                     "Al-Ketan O., Abu Al-Rub R.K. Multifunctional mechanical metamaterials based on TPMS. Adv. Eng. Mater. 21 (2019)",
                 status="по памяти; закон жёсткости и прочности гироида — свой расчёт organic/structure/lattice_results.json"),
    "voronoi_lattice": dict(ref="Gibson L.J., Ashby M.F. Cellular Solids. Cambridge (1997) — пены и решётки; открытая ячейка "
                                "Вороного — общий приём печатных решёток", status="по памяти"),
    "murray": dict(ref="Murray C.D. The physiological principle of minimum work. PNAS 12:207–214 (1926)",
                   status="по памяти; показатель 2,41 для турбулентного газа — свой расчёт organic/structure/manifold_results.json"),
    "lsystem": dict(ref="Lindenmayer A. Mathematical models for cellular interactions in development. J. Theor. Biol. 18 (1968)",
                    status="по памяти"),
    "tubercles": dict(ref="Fish F.E., Battle J.M. Hydrodynamic design of the humpback whale flipper. J. Morphol. 225:51–60 (1995); "
                          "Miklosovic D.S. et al. Leading-edge tubercles delay stall on humpback whale flippers. Phys. Fluids 16:L39 (2004); "
                          "Johari H. et al. Effects of leading-edge protuberances on airfoil performance. AIAA J. 45 (2007)",
                      status="по памяти; задержка срыва — для изолированного крыла при малых числах Маха; для решётки "
                             "лопаток ГТУ — гипотеза"),
    "serrations": dict(ref="Howe M.S. Noise produced by a sawtooth trailing edge. J. Acoust. Soc. Am. 90(1):482–487 (1991)",
                       status="по памяти; шум — да, КПД ступени — гипотеза"),
    "riblets": dict(ref="Walsh M.J. Riblets as a viscous drag reduction technique. AIAA J. 21(4) (1983); Bechert D.W. et al. "
                        "Experiments on drag-reducing surfaces and their optimization with an adjustable geometry. J. Fluid Mech. 338 (1997)",
                    status="по памяти; снижение трения до ~5–8 % при шаге s⁺ ≈ 15–20 — порядок по памяти, оценка"),
    "bow_lean_sweep": dict(ref="навал, изгиб и стреловидность лопаток — Denton J.D., Xu L. The exploitation of 3D flow in turbomachinery "
                               "design. VKI LS 1999-02; Gümmer V. et al. Using sweep and dihedral to control three-dimensional flow "
                               "in transonic stators. J. Turbomach. 123 (2001)", status="по памяти"),
    "shell_vibration": dict(ref="Leissa A.W. Vibration of Shells. NASA SP-288 (1973) — уравнения Доннелла; Soedel W. Vibrations of "
                                "Shells and Plates. Dekker (2004)", status="по памяти; формула проверена тестом на изотропном пределе "
                                "и расчётом CalculiX (семейство shell)"),
    "smeared_stiffeners": dict(ref="Baruch M., Singer J. Effect of eccentricity of stiffeners on the general instability of stiffened "
                                   "cylindrical shells under hydrostatic pressure. J. Mech. Eng. Sci. 5(1):23–27 (1963)", status="по памяти"),
    "victoria_ribs": dict(ref="рёбра листа виктории амазонской (Victoria amazonica) как прообраз ребристых конструкций (оранжерея и "
                              "Хрустальный дворец Дж. Пакстона, 1850-е) — общеизвестный пример биомиметики", status="по памяти; "
                              "выигрыш жёсткости ветвящихся рёбер — гипотеза, меряется расчётом"),
}
