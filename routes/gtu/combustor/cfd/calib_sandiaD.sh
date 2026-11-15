#!/bin/bash
# Эталон Sandia Flame D (RANS) до уровня 1 воронки камеры: та же версия OpenFOAM и те же модели, что кейсы L1
# (reactingFoam LTS, k-epsilon, EDC v2005, P1, полный GRI-Mech 3.0 + TDAC) - навык openfoam-combustor, разделы 5-6.
# Запуск из корня клона в окружении OpenFOAM (контейнер opencfd/openfoam-default:2512):
#   bash routes/gtu/combustor/cfd/calib_sandiaD.sh 1.44     # C1 по умолчанию; второй прогон - 1.6 (чувствительность)
# Пишет в $GTU_OUT/SandiaD_C1_<C1>/: кейс, журналы, профили sample (оси и станции z/D 7,5 ... 75).
set -euo pipefail
C1="${1:-1.44}"
OUT="${GTU_OUT:-/tmp/sandia}"
D="$OUT/SandiaD_C1_$C1"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$D"
T="${FOAM_TUTORIALS:-/nonexistent}/combustion/reactingFoam/RAS/SandiaD_LTS"
if [ -d "$T" ]; then
  cp -r "$T/." "$D/"
else    # в образе нет учебников - те же файлы v2512 из исходников (gitlab.com/openfoam/core/openfoam, тег OpenFOAM-v2512)
  R=https://gitlab.com/openfoam/core/openfoam/-/raw/OpenFOAM-v2512/tutorials/combustion/reactingFoam/RAS/SandiaD_LTS
  for f in Allrun 0.orig/CH4 0.orig/CO 0.orig/CO2 0.orig/G 0.orig/H 0.orig/H2 0.orig/H2O 0.orig/N2 0.orig/O 0.orig/O2 \
           0.orig/OH 0.orig/T 0.orig/U 0.orig/Ydefault 0.orig/alphat 0.orig/epsilon 0.orig/k 0.orig/nut 0.orig/p \
           constant/boundaryRadiationProperties constant/chemistryProperties constant/combustionProperties \
           constant/fvOptions constant/g constant/radiationProperties constant/thermophysicalProperties \
           constant/turbulenceProperties system/blockMeshDict system/controlDict system/fvSchemes system/fvSolution \
           system/sampleDict system/setFieldsDict; do
    mkdir -p "$D/$(dirname "$f")"; curl -sSf "$R/$f" -o "$D/$f"
  done
fi
# полный GRI-Mech 3.0 (53/325, с азотом) вместо урезанного grimech30.dat учебника (36/219, без NO) - навык, разделы 5, 9
cp "$HERE/foam/reactionsGRI" "$HERE/foam/thermo.compressibleGasGRI" "$D/constant/"
sed -i '/runApplication chemkinToFoam/,+2d' "$D/Allrun"
# опечатка образца v2606/v2512 (набор Radial_15: две компоненты) и NO, CO в выборке - навык, разделы 6, 9
sed -i 's/end (0.024 0.108);/end (0.024 0 0.108);/' "$D/system/sampleDict"
sed -i 's/^\(\s*fields\s*(\)/\1 NO CO /' "$D/system/sampleDict"
if [ "$C1" != "1.44" ]; then
  foamDictionary "$D/constant/turbulenceProperties" -entry RAS.kEpsilonCoeffs.C1 -set "$C1" > /dev/null
fi
chmod +x "$D/Allrun"
( cd "$D" && ./Allrun && postProcess -dict system/sampleDict -latestTime > log.sample 2>&1 )
tail -n 3 "$D"/log.* > "$D/_tails.txt" || true
