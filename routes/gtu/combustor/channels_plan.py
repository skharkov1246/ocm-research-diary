"""Каналы камеры KV и печатные внутренние структуры: постановка сопряжённой (adjoint) оптимизации формы и
топологической оптимизации, число параметров, цена итерации, смета -> channels_plan_results.json и шаблоны в cfd/.

Каналы (README, раздел «Каналы»):
  * диффузор за второй ступенью компрессора до кожуха камеры - SU2 v8.6.0, сжимаемый RANS SST, осесимметрия
    (AXISYMMETRIC= YES), дискретное сопряжённое (MATH_PROBLEM= DISCRETE_ADJOINT), FFD 2D (FFD_CONTROL_POINT_2D);
    цель - потери полного давления + неравномерность на входе в камеру (SURFACE_TOTAL_PRESSURE, SURFACE_UNIFORMITY);
  * переходный канал жаровая труба -> кольцо соплового аппарата турбины и улитка (сбор воздуха к выносной камере) -
    OpenFOAM openfoam.com adjointOptimisationFoam (несжимаемый: число Маха 0,1-0,25 - тренды формы; итог проверяется
    сжимаемым SU2), цели PtLosses + uniformityPatch, объёмные B-сплайны (volumetricBSplines), BFGS.
Печатные структуры: топливные коллекторы (равная раздача на n элементов) - adjointOptimisationFoam, топологические
переменные (topODesignVariables), цель flowRatePartition + PtLosses; решётки охлаждения - параметрическая
ячейка TPMS по CHT-расчёту ячейки (chtMultiRegionFoam) и модели-заменителю, топологическая постановка - pyMOTO
(теплопроводность и термомеханика); крепления картриджа - TopOpt_in_PETSc или SU2 (TOPOLOGY_OPTIMIZATION).
Все имена ключей и объектов сверены по исходникам (SU2 v8.6.0 CConfig.cpp и option_structure.hpp; openfoam.com
тег OpenFOAM-v2512: дерево src/optimisation/adjointOptimisation, учебники sbend и fork-uneven) - раздел sources.

Цены итерации: SU2 - 1,1e-5 с на узел·итерацию на ядро (из сметы пакета SU2 KV-3: s_per_iter_base 0,782 с при
1 137 992 узлах на 16 процессах, platform/su2/jobs_aws_su2_kv3.json); OpenFOAM simpleFoam k-omega SST - замер этой
сессии (apt-1912, 24 000 ячеек, 300 итераций, 4 общих ядра под нагрузкой); сопряжённое - 1,5 прямого (оценка).
Запуск: python channels_plan.py [--cconfig /path/CConfig.cpp]  (без файла ключи не сверяются заново)
"""
import argparse, json, math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
DATE = "2026-11-14"
PRICE = dict(expected=0.7064, planning=0.75)       # долл./ч, как сметы уровней (README, «Смета»)
VCPU = 64
SU2_S_NODE_IT = 0.782 * 16 / 1137992               # с на узел·итерацию на ядро
SF_S_CELL_IT = (9.05 - 0.19) / 299 / 24000          # замер simpleFoam этой сессии
ADJ_FACTOR = 1.5                                    # сопряжённое / прямое (оценка)

SU2_TPL = """% Диффузор за ступенью 2 компрессора KV-3 -> кожух камеры: сопряжённая оптимизация формы (SU2 v8.6.0)
% Шаблон: <...> заполняет генератор сетки канала (не сделан - README, «Каналы»). Ключи сверены с CConfig.cpp v8.6.0.
SOLVER= RANS
KIND_TURB_MODEL= SST
MATH_PROBLEM= DISCRETE_ADJOINT
AXISYMMETRIC= YES
RESTART_SOL= NO
MACH_NUMBER= 0.3
INIT_OPTION= TD_CONDITIONS
FREESTREAM_PRESSURE= <p_stat_inlet>
FREESTREAM_TEMPERATURE= <T_stat_inlet>
REF_DIMENSIONALIZATION= DIMENSIONAL
MARKER_INLET= ( inlet, <T0_inlet>, <p0_inlet>, <dir_r>, <dir_t>, 0.0 )
MARKER_OUTLET= ( outlet, <p_outlet> )
MARKER_HEATFLUX= ( hub, 0.0, shroud, 0.0 )
MARKER_SYM= ( axis )
MARKER_ANALYZE= ( outlet, inlet )
MARKER_ANALYZE_AVERAGE= MASSFLUX
OBJECTIVE_FUNCTION= SURFACE_TOTAL_PRESSURE, SURFACE_UNIFORMITY
OBJECTIVE_WEIGHT= 1.0, <w_uniformity>
DV_KIND= FFD_CONTROL_POINT_2D
DV_MARKER= ( hub, shroud )
DV_PARAM= ( BOX, 0, 0, 0.0, 1.0 )
FFD_DEFINITION= ( BOX, <x0>, <y0>, 0.0, <x1>, <y0>, 0.0, <x1>, <y1>, 0.0, <x0>, <y1>, 0.0, 0,0,0, 0,0,0, 0,0,0, 0,0,0 )
FFD_DEGREE= ( <nx>, <ny>, 0 )
FFD_CONTINUITY= 2ND_DERIVATIVE
DEFORM_MESH= YES
DEFORM_LINEAR_SOLVER= FGMRES
DEFORM_LINEAR_SOLVER_ITER= 500
DEFORM_STIFFNESS_TYPE= WALL_DISTANCE
OPT_OBJECTIVE= SURFACE_TOTAL_PRESSURE * 1.0
OPT_CONSTRAINT= NONE
OPT_ITERATIONS= <n_design>
OPT_ACCURACY= 1E-6
OPT_BOUND_UPPER= 0.01
OPT_BOUND_LOWER= -0.01
OPT_GRADIENT_FACTOR= 1.0
ITER= 4000
CONV_FIELD= RMS_DENSITY
CONV_RESIDUAL_MINVAL= -10
MESH_FILENAME= <mesh.su2>
MESH_FORMAT= SU2
"""

OF_TPL = """// Переходный канал (или улитка) KV: adjointOptimisationFoam openfoam.com v2512 - шаблон optimisationDict.
// Объекты и ключи - как в учебниках v2512 sbend (PtLosses, volumetricBSplines, BFGS) и fork-uneven (uniformityPatch).
optimisationManager     steadyOptimisation;
primalSolvers
{
    p1 { active true; type incompressible; solver simple;
         solutionControls { nIters 3000; residualControl { "p.*" 1e-7; "U.*" 1e-7; } } }
}
adjointManagers
{
    am1
    {
        primalSolver p1;
        adjointSolvers
        {
            as1
            {
                active true; type incompressible; solver adjointSimple;
                objectives
                {
                    type incompressible;
                    objectiveNames
                    {
                        losses     { weight 1;   type PtLosses;        patches (inlet outlet); }
                        uniformity { weight <w>; type uniformityPatch; patches (outlet); }
                    }
                }
                ATCModel { ATCModel standard; }
                solutionControls { nIters 3000; residualControl { "pa.*" 1e-7; "Ua.*" 1e-7; } }
            }
        }
    }
}
optimisation
{
    designVariables { type shape; shapeType volumetricBSplines; sensitivityType shapeFI; patches (walls); maxInitChange 2e-3; }
    updateMethod { method BFGS; BFGS { etaHessian 0.8; nSteepestDescent 1; scaleFirstHessian true; } }
}
// constant/dynamicMeshDict: solver volumetricBSplinesMotionSolver; box: type cylindrical; nCPsU <nu>; nCPsV <nv>; nCPsW <nw>
"""

SU2_KEYS_USED = sorted(set(re.findall(r"^([A-Z][A-Z0-9_]+)=", SU2_TPL, flags=re.M)))
SU2_ENUM_VALUES = ["SURFACE_TOTAL_PRESSURE", "SURFACE_UNIFORMITY", "FFD_CONTROL_POINT_2D", "DISCRETE_ADJOINT", "RANS",
                   "SST", "MASSFLUX", "TOPOL_COMPLIANCE", "VOLUME_FRACTION", "HEAT_EQUATION", "MULTIPHYSICS", "INC_RANS"]
OF_VERIFIED = dict(
    solver="applications/solvers/incompressible/adjointOptimisationFoam (adjointOptimisationFoam.C)",
    objectives_incompressible=["objectiveFlowRate", "objectiveFlowRatePartition", "objectiveForce", "objectiveMoment",
                               "objectiveNutSqr", "objectivePowerDissipation", "objectivePtLosses",
                               "objectiveUniformityCellZone", "objectiveUniformityPatch"],
    design_variables=["designVariables", "levelSet", "shape", "topODesignVariables"],
    parameterization=["Bezier", "NURBS (NURBS3DVolume: cartesian, cylindrical, controlPointsDefinition, volBSplinesBase)"],
    tutorials=["shapeOptimisation/sbend (PtLosses, volumetricBSplines transformBox nCPs 9×5×3, BFGS, 3000 + 3000 итераций)",
               "shapeOptimisation/fork-uneven (flowRate, flowRatePartition, uniformityPatch)",
               "topologyOptimisation/monoFluidAero"],
    note="тепловой цели (температура, тепловой поток) в несжимаемых объектах v2512 нет - сопряжённого теплообмена "
         "в adjointOptimisationFoam нет; для решёток охлаждения - другой путь (ниже)")


def check_su2_keys(path):
    if not path or not os.path.exists(path):
        return None
    t = open(path, encoding="utf-8", errors="replace").read()
    keys = {k: bool(re.search(r'add[A-Za-z]+Option\("%s"' % k, t)) for k in SU2_KEYS_USED}
    return dict(file=path, all_present=all(keys.values()), keys=keys)


def check_su2_enums(path, cconfig=None):
    """значения перечислений: MakePair в option_structure.hpp; MATH_PROBLEM разбирается своим классом опции -
    DISCRETE_ADJOINT сверяется по сообщению CConfig.cpp («MATH_PROBLEM= DISCRETE_ADJOINT requires AD support»)"""
    if not path or not os.path.exists(path):
        return None
    t = open(path, encoding="utf-8", errors="replace").read()
    out = {v: f'MakePair("{v}"' in t for v in SU2_ENUM_VALUES}
    if cconfig and os.path.exists(cconfig):
        out["DISCRETE_ADJOINT"] = "MATH_PROBLEM= DISCRETE_ADJOINT" in open(cconfig, encoding="utf-8", errors="replace").read()
    return out


def cost(core_h, label):
    mh = core_h / VCPU
    return dict(core_hours=core_h, machine_hours=mh, usd_expected=mh * PRICE["expected"], usd_planning=1.5 * mh * PRICE["planning"],
                basis=label)


def plan():
    f3 = json.load(open(os.path.join(GTU, "frame3", "frame3_point_results.json"), encoding="utf-8"))
    st2 = f3["stage2"]
    tz = f3["turbine"]["options"]["3"]
    D2 = 0.703        # frame3 README, таблица ступеней: D2 ступени 2 = 0,703 м (frame3_geom: stages.stage2)
    b2 = st2["b2_m"]
    Dm, h1 = tz["Dm"], tz["stages"][0]["h_rotor_in_m"]
    ch = []
    # диффузор: 2D осесимметрия, FFD 12 × 4 точек, по 2 направления; закреплены ряды на входе и выходе
    n_ffd, fixed = 12 * 4 * 2, 2 * 4 * 2
    nodes = 60000
    it_primal, it_adj, n_design = 4000, 4000, 30
    ci = nodes * (it_primal + ADJ_FACTOR * it_adj) * SU2_S_NODE_IT / 3600
    ch.append(dict(id="diffuser", what="диффузор за ступенью 2 (D2 = 0,703 м, b2 = %.1f мм) до кожуха камеры" % (1e3 * b2),
                   tool="SU2 v8.6.0, RANS SST, AXISYMMETRIC, DISCRETE_ADJOINT, FFD_CONTROL_POINT_2D",
                   objective="потери полного давления (SURFACE_TOTAL_PRESSURE) + w·неравномерность (SURFACE_UNIFORMITY) на выходе",
                   build_note="дискретное сопряжённое требует сборки SU2 с автодифференцированием (CConfig.cpp: «MATH_PROBLEM= "
                              "DISCRETE_ADJOINT requires AD support»); install_su2.sh узла su2 собирает без него - "
                              "нужна отдельная сборка с -Denable-autodiff=true (meson_options.txt v8.6.0: option enable-autodiff, по умолчанию false - открыт)",
                   constraints=["площадь выхода = площади входа в кожух (не меньше)", "длина не больше габарита frame3",
                                "угол стенки ≤ 8° на участке торможения (оценка)", "минимальная толщина стенки печати 1 мм"],
                   design_variables=dict(ffd_points="12 × 4", directions=2, total=n_ffd, active=n_ffd - fixed),
                   mesh_nodes=nodes, iterations=dict(primal=it_primal, adjoint=it_adj, design=n_design),
                   core_hours_per_design_iter=ci, **cost(ci * n_design, "SU2 1,1e-5 с/(узел·итер.·ядро) × (прямое + 1,5·сопр.)"),
                   inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json", field="stage2.b2_m", used=b2)]))
    # переходный канал: 3D сектор (одна труба или сектор кольца), B-сплайны 9 × 5 × 3 (как учебник sbend), 3 направления
    for cid, what, cells, cps, n_design in (
            ("transition", "переходный канал жаровая труба -> кольцо соплового аппарата (Dm = %.3f м, h = %.1f мм)"
             % (Dm, 1e3 * h1), 1.5e6, (9, 5, 3), 40),
            ("scroll", "улитка сбора воздуха за диффузором к выносной камере (вход картриджа)", 3.0e6, (12, 6, 6), 40)):
        tot = cps[0] * cps[1] * cps[2] * 3
        active = (cps[0] - 2) * cps[1] * cps[2] * 3        # закреплены два крайних ряда по оси (стыки)
        ci = cells * 3000 * (1 + ADJ_FACTOR) * SF_S_CELL_IT / 3600
        ch.append(dict(id=cid, what=what, tool="openfoam.com adjointOptimisationFoam (v2512), несжимаемый, k-omega SST "
                                               "(сопряжённая турбулентность - замороженная, ATC standard)",
                       objective="PtLosses (вход-выход) + w·uniformityPatch (выход)",
                       constraints=["площади и положения стыков закреплены (крайние ряды точек)",
                                    "минимальный радиус кривизны и толщина стенки печати", "без изменения длины модуля"],
                       design_variables=dict(control_points=" × ".join(map(str, cps)), directions=3, total=tot,
                                             active=active),
                       mesh_cells=cells, iterations=dict(primal=3000, adjoint=3000, design=n_design),
                       core_hours_per_design_iter=ci,
                       **cost(ci * n_design, "simpleFoam 1,23e-6 с/(ячейка·итер.·ядро) - замер; сопряжённое 1,5 прямого (оценка)"),
                       inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json",
                                    field="turbine.options.3.Dm, stages[0].h_rotor_in_m", used=[Dm, h1])] if cid == "transition"
                       else []))
    topo = []
    cells = 1.0e6
    ci = cells * 3000 * (1 + ADJ_FACTOR) * SF_S_CELL_IT / 3600
    topo.append(dict(id="fuel_manifold", what="печатный топливный коллектор: один вход (один клапан) -> n элементов, "
                                              "равная раздача на 30-100 % расхода",
                     tool="adjointOptimisationFoam v2512, topODesignVariables (пористость Бринкмана по ячейкам)",
                     objective="flowRatePartition (заданные доли на выходах) + w·PtLosses",
                     design_variables=dict(kind="значение на ячейку", total=int(cells)), cycles=150,
                     core_hours_per_cycle=ci, **cost(ci * 150, "simpleFoam замер × (1 + 1,5) × 150 циклов (оценка)"),
                     manufacturability="фильтр радиуса ≥ 2 ячеек и порог минимального канала 0,3 мм (comb_space.D_JET_MIN, "
                                       "печать); нависание ≤ 45° и удаление порошка - проверка после (не в постановке)"))
    pts, cells_u = 40, 2.0e5
    ci = cells_u * 3000 * 1.2 * 6.0e-5 / 3600 / 1.0
    topo.append(dict(id="cooling_lattice", what="печатная решётка охлаждения стенки картриджа (двойная стенка)",
                     tool="этап 1 - параметрическая ячейка TPMS (гироид, алмаз: размер ячейки, толщина, пористость) - "
                          "chtMultiRegionFoam ячейки с периодикой + модель-заменитель (platform/surrogate); этап 2 - "
                          "топология проводимости и термомеханики - pyMOTO (открыт, README: «Thermal and thermo-mechanic "
                          "coupling»), сопряжённый теплообмен течения - пакета с готовой постановкой не найдено",
                     objective="наименьшая температура металла при перепаде dp_wall ≤ 1 % p3 и расходе охлаждения a_cool",
                     design_variables=dict(stage1="3-4 параметра ячейки", stage2="плотность на элемент"),
                     points=pts, core_hours_per_point=ci,
                     **cost(ci * pts, "CHT ячейки 0,2 млн × 3000 итер. × 1,2·6,0e-5 с (реагирующий замер L1 как верхняя "
                                      "граница для нереагирующего) × 40 точек (оценка)")))
    topo.append(dict(id="cartridge_bracket", what="крепления и фланец картриджа (жёсткость, масса)",
                     tool="TopOpt_in_PETSc (DTU, открыт, README) или SU2 TOPOLOGY_OPTIMIZATION (TOPOL_COMPLIANCE, "
                          "VOLUME_FRACTION, ключи v8.6.0), проверка - CalculiX (навык calculix-rotor)",
                     objective="наименьшая податливость при доле объёма; ограничение - частоты вне полосы пульсаций",
                     **cost(20.0, "TopOpt на 1-2 млн элементов: десятки ядро-часов (оценка)")))
    return ch, topo


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cconfig", help="CConfig.cpp SU2 v8.6.0 для сверки ключей")
    ap.add_argument("--options", help="option_structure.hpp SU2 v8.6.0 для сверки значений")
    a = ap.parse_args(argv)
    os.makedirs(os.path.join(HERE, "cfd"), exist_ok=True)
    open(os.path.join(HERE, "cfd", "diffuser_adjoint.cfg.tpl"), "w", encoding="utf-8").write(SU2_TPL)
    open(os.path.join(HERE, "cfd", "channel_optimisationDict.tpl"), "w", encoding="utf-8").write(OF_TPL)
    ch, topo = plan()
    out = os.path.join(HERE, "channels_plan_results.json")
    prev = json.load(open(out, encoding="utf-8")) if os.path.exists(out) else {}
    su2 = check_su2_keys(a.cconfig) or prev.get("su2_key_check")
    en = check_su2_enums(a.options, a.cconfig) or prev.get("su2_enum_check")
    res = dict(title="Каналы и печатные структуры камеры KV-3: постановка сопряжённой и топологической оптимизации",
               method="постановка и смета; цены итерации - из замеров и смет (см. basis в каждой строке)",
               software="python; сверка ключей - исходники SU2 v8.6.0 и openfoam.com v2512", date=DATE,
               status="план (расчётов оптимизации не было); числа сметы - оценка",
               throughput=dict(su2_s_per_node_iter_core=SU2_S_NODE_IT,
                               su2_basis="routes/gtu/platform/su2/jobs_aws_su2_kv3.json → estimate.jobs[0].s_per_iter_base "
                                         "0,782 с при 1 137 992 узлах, 16 процессов (смета пакета SU2, не замер этой сессии)",
                               simplefoam_s_per_cell_iter_core=SF_S_CELL_IT,
                               simplefoam_basis="замер 14.11.2026: apt OpenFOAM 1912, канал 400×60 ячеек, k-omega SST, 300 "
                                                "итераций, ExecutionTime 9,05 − 0,19 с; 4 общих ядра под нагрузкой 23-28",
                               adjoint_factor=ADJ_FACTOR, adjoint_factor_basis="оценка"),
               channels=ch, topology=topo, su2_keys_used=SU2_KEYS_USED, su2_key_check=su2, su2_enum_check=en,
               openfoam_adjoint_verified=OF_VERIFIED,
               totals=dict(core_hours=sum(x["core_hours"] for x in ch + topo),
                           usd_expected=sum(x["usd_expected"] for x in ch + topo),
                           usd_planning=sum(x["usd_planning"] for x in ch + topo)),
               sources=["https://raw.githubusercontent.com/su2code/SU2/v8.6.0/Common/src/CConfig.cpp (открыт)",
                        "https://raw.githubusercontent.com/su2code/SU2/v8.6.0/Common/include/option_structure.hpp (открыт)",
                        "https://gitlab.com/openfoam/core/openfoam (тег OpenFOAM-v2512; дерево adjointOptimisation и "
                        "учебники - открыт через API GitLab)",
                        "https://raw.githubusercontent.com/topopt/TopOpt_in_PETSc/master/README.md (открыт)",
                        "https://raw.githubusercontent.com/aatmdelissen/pyMOTO/master/README.md (открыт)",
                        "https://raw.githubusercontent.com/dolfin-adjoint/pyadjoint/master/README.md (открыт)"])
    json.dump(res, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(res["totals"], ensure_ascii=False), "ключи SU2:", (su2 or {}).get("all_present"), "значения:", en)


if __name__ == "__main__":
    main()
