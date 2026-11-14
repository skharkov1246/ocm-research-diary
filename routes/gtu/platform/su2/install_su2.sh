#!/usr/bin/env bash
# Сборка SU2 v8.6.0 «Harrier» с MPI из исходников (пакетов SU2 в apt Ubuntu 24.04 и на PyPI нет — см. README).
# Проверено 14.11.2026: Ubuntu 24.04, g++ 13.3, OpenMPI 4.1.6, 4 ядра — setup ~2 мин, сборка ~8 мин.
# Использование:  bash install_su2.sh [префикс установки] [папка исходников]
#   по умолчанию префикс $HOME/su2-8.6.0-mpi, исходники $HOME/SU2-src.
# После установки:  export SU2_RUN=<префикс>/bin PATH=<префикс>/bin:$PATH PYTHONPATH=<префикс>/bin:$PYTHONPATH
# meson.py сам скачивает meson, ninja и внешние библиотеки (подмодули) с GitHub — нужен доступ к github.com.
set -euo pipefail
PREFIX="${1:-$HOME/su2-8.6.0-mpi}"
SRC="${2:-$HOME/SU2-src}"
JOBS="${JOBS:-$(nproc)}"

if ! command -v mpicxx >/dev/null; then
  SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo"
  $SUDO apt-get update -y && $SUDO apt-get install -y --no-install-recommends \
    build-essential python3 git pkg-config openmpi-bin libopenmpi-dev
fi
if [ ! -d "$SRC/.git" ]; then
  git clone --depth 1 -b v8.6.0 https://github.com/su2code/SU2.git "$SRC"
fi
cd "$SRC"
if [ -d build_mpi ]; then   # повторный запуск: перенастроить существующую папку сборки
  ./meson.py setup --reconfigure build_mpi -Dwith-mpi=enabled -Dbuildtype=release --prefix="$PREFIX"
else
  ./meson.py setup build_mpi -Dwith-mpi=enabled -Dbuildtype=release --prefix="$PREFIX"
fi
./ninja -C build_mpi -j"$JOBS" install
"$PREFIX/bin/SU2_CFD" --help | head -1
echo "SU2_RUN=$PREFIX/bin"
