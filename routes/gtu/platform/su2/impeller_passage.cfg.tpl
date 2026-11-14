% SU2 v8.6.0: один межлопаточный канал центробежного колеса (основная лопатка + полулопатка),
% сетка routes/gtu/platform/geom/impeller.py. Файл заполняет run_case.py: поля $$ИМЯ -> числа из JSON сетки
% и ключей командной строки. Имена ключей сверяются с Common/src/CConfig.cpp v8.6.0 (run_case.py --check-keys).
% Схема и пределы — см. README.md этой папки; что проверено прогоном, а что нет — там же.
%
% ---- задача и модель ----
SOLVER= RANS
KIND_TURB_MODEL= ${TURB}
${SA_OPTIONS}
MATH_PROBLEM= DIRECT
RESTART_SOL= ${RESTART}
SOLUTION_FILENAME= ${SOLUTION_FILENAME}
%
% ---- начальное поле и набегающий поток (размерный расчёт) ----
% начальная скорость: MACH_NUMBER по направлению AOA (AOA= 90 -> вдоль +z, оси машины; код CEulerSolver.cpp)
MACH_NUMBER= ${INIT_MACH}
AOA= ${AOA}
SIDESLIP_ANGLE= 0.0
FREESTREAM_OPTION= TEMPERATURE_FS
FREESTREAM_PRESSURE= ${P01}
FREESTREAM_TEMPERATURE= ${T01}
FREESTREAM_TURBULENCEINTENSITY= ${TU}
FREESTREAM_TURB2LAMVISCRATIO= ${MUT_RATIO}
INIT_OPTION= TD_CONDITIONS
REF_DIMENSIONALIZATION= DIMENSIONAL
%
% ---- газ: воздух, идеальный газ ----
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= ${GAMMA}
GAS_CONSTANT= ${R_GAS}
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716E-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.4
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
%
% ---- вращающаяся система координат (ось z, рад/с; знак — из JSON сетки) ----
GRID_MOVEMENT= ROTATING_FRAME
MOTION_ORIGIN= 0.0 0.0 0.0
ROTATION_RATE= ${OMEGA_X} ${OMEGA_Y} ${OMEGA_Z}
MACH_MOTION= 0.35
% рампа частоты: (начальная частота рад/с, обновлять каждые N итераций, длина рампы в итерациях);
% рампа меняет только z-компоненту (код CConfig::SetRate) и работает только с MARKER_TURBOMACHINERY
RAMP_MOTION_FRAME= ${RAMP_ROT}
RAMP_MOTION_FRAME_COEFF= (${RAMP_ROT_INI}, ${RAMP_ROT_FREQ}, ${RAMP_ROT_ITER})
%
% ---- стенки: адиабатные; покрышка неподвижна в абсолютной системе (MARKER_SHROUD) ----
MARKER_HEATFLUX= (hub, 0.0, blade, 0.0, shroud, 0.0)
MARKER_SHROUD= (shroud)
%
% ---- периодичность поворотом: periodic_2 = periodic_1, повёрнутый на +360/Z вокруг +z ----
MARKER_PERIODIC= (${PER_MASTER}, ${PER_SLAVE}, 0.0, 0.0, 0.0, 0.0, 0.0, ${ANGLE}, 0.0, 0.0, 0.0)
%
% ---- вход и выход: условия Giles на турбомашинных маркерах ----
MARKER_TURBOMACHINERY= (inlet, outlet)
TURBOMACHINERY_KIND= AXIAL_CENTRIFUGAL
TURBO_PERF_KIND= COMPRESSOR
% Giles: маркер, тип, p0 [Па], T0 [К], направление в осях решётки (нормаль, окружное, по высоте),
% релаксация среднего, релаксация Фурье; выход: статическое давление [Па]
MARKER_GILES= (inlet, TOTAL_CONDITIONS_PT, ${P01}, ${T01}, 1.0, 0.0, 0.0, ${RELAX_IN}, 0.0, outlet, STATIC_PRESSURE_1D, ${P_OUT}, 0.0, 0.0, 0.0, 0.0, ${RELAX_OUT}, 0.0)
% рампа давления на выходе (Giles): (начальное давление Па, обновлять каждые N итераций, длина рампы)
RAMP_OUTLET= ${RAMP_OUT}
RAMP_OUTLET_COEFF= (${RAMP_OUT_INI}, ${RAMP_OUT_FREQ}, ${RAMP_OUT_ITER})
AVERAGE_PROCESS_KIND= MIXEDOUT
PERFORMANCE_AVERAGE_PROCESS_KIND= MIXEDOUT
SPANWISE_KIND= EQUISPACED
NUM_SPANWISE_SECTIONS= ${NSPAN}
MARKER_ANALYZE= (inlet, outlet)
MARKER_ANALYZE_AVERAGE= MASSFLUX
%
% ---- численная схема ----
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
CONV_NUM_METHOD_FLOW= ${FLOW_SCHEME}
MUSCL_FLOW= ${MUSCL}
SLOPE_LIMITER_FLOW= ${LIMITER}
VENKAT_LIMITER_COEFF= ${VENKAT}
ENTROPY_FIX_COEFF= ${ENTROPY_FIX}
JST_SENSOR_COEFF= (0.5, 0.02)
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_NUMBER= ${CFL}
CFL_ADAPT= ${CFL_ADAPT}
CFL_ADAPT_PARAM= (${CFL_ADAPT_PARAM})
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ${LIN_PREC}
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= ${LIN_ITER}
MGLEVEL= 0
%
% ---- останов ----
ITER= ${ITER}
CONV_FIELD= RMS_DENSITY
CONV_RESIDUAL_MINVAL= ${RES_MIN}
CONV_STARTITER= 10
%
% ---- вывод ----
SCREEN_OUTPUT= (INNER_ITER, WALL_TIME, RMS_DENSITY, RMS_ENERGY, RMS_TKE, SURFACE_MASSFLOW, SURFACE_TOTAL_PRESSURE)
HISTORY_OUTPUT= (ITER, WALL_TIME, RMS_RES, FLOW_COEFF, FLOW_COEFF_SURF, TURBO_PERF)
SCREEN_WRT_FREQ_INNER= 10
VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)
OUTPUT_FILES= (${OUTPUT_FILES})
OUTPUT_WRT_FREQ= ${WRT_FREQ}
TABULAR_FORMAT= CSV
CONV_FILENAME= ${HISTORY_NAME}
RESTART_FILENAME= ${RESTART_NAME}
VOLUME_FILENAME= ${VOLUME_NAME}
MESH_FILENAME= ${MESH}
MESH_FORMAT= SU2
