# syntax=docker/dockerfile:1.7
# Образ платформы трека «Простая литая ГТУ-заменитель» (issue #59) для AWS Batch и локального запуска:
# Python 3.12 + cantera, numpy, scipy, gmsh, boto3; CalculiX ccx 2.21; SU2 v8.6.0 «Harrier» с MPI; код трека.
# Точка входа — routes/gtu/platform/batch/entrypoint.py: запуск задания по JSON (README узла batch).
#
# Сборка — из корня ЧИСТОЙ копии коммита (git worktree add --detach), чтобы в образ не попали чужие
# незакоммиченные файлы, и с хешем этого коммита:
#   docker build -f routes/gtu/platform/docker/Dockerfile \
#       --build-arg GIT_SHA="$(git rev-parse HEAD)" -t gtu-platform:"$(git rev-parse --short=12 HEAD)" .
#
# Почему SU2 собирается, а не берётся готовым: у выпуска v8.6.0 нет архивов под именами прошлых выпусков
# (SU2-v8.6.0-linux64.zip, -linux64-mpi.zip, -linux64-omp.zip — 404; у v8.5.0 есть -linux64-mpi и -linux64-omp,
# проверено 14.11.2026 запросами к github.com/su2code/SU2/releases/download), а сборка из исходника с MPI
# уже проверена в узле su2 (install_su2.sh: Ubuntu 24.04, g++ 13.3, OpenMPI 4.1.6). Своя сборка к тому же
# связывает SU2 с той же OpenMPI, что стоит в образе. Коммит тега сверяется при сборке.
#
# Основа — Ubuntu 24.04, а не python:3.12-slim (Debian 12): системный Python 3.12, calculix-ccx 2.21 (как в
# навыке calculix-rotor) и OpenMPI 4.1.6 — те же, на которых проверены узлы платформы.
# BASE_IMAGE можно подменить зеркалом (Docker Hub ограничивает анонимные скачивания: 429).

ARG BASE_IMAGE=ubuntu:24.04

# ---------- стадия 1: сборка SU2 ----------
FROM ${BASE_IMAGE} AS su2-build
ARG SU2_TAG=v8.6.0
ARG SU2_COMMIT=7048fd35f56d4bf52e586429d6746a36a10255f8
ARG BUILD_JOBS=4
ENV DEBIAN_FRONTEND=noninteractive
# Версии apt не закрепляются (DL3008): архив Ubuntu убирает старые сборки после обновлений безопасности, и
# закреплённая сборка перестанет собираться; фактические версии пишутся в /opt/gtu/PACKAGES.txt образа.
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential python3 git pkg-config ca-certificates openmpi-bin libopenmpi-dev \
 && rm -rf /var/lib/apt/lists/*
# extra_ca — необязательный секрет сборки: набор корневых сертификатов для сети с TLS-прокси.
# Действует только внутри этой команды и в слой образа не попадает; без него — системные сертификаты.
# cd внутри одной команды (DL3003): клон, setup и сборке нужен один и тот же необязательный секрет.
# hadolint ignore=DL3003
RUN --mount=type=secret,id=extra_ca,required=false \
    set -eu; \
    if [ -s /run/secrets/extra_ca ]; then \
      export GIT_SSL_CAINFO=/run/secrets/extra_ca SSL_CERT_FILE=/run/secrets/extra_ca; fi; \
    git clone --depth 1 -b "$SU2_TAG" https://github.com/su2code/SU2.git /src/SU2; \
    test "$(git -C /src/SU2 rev-parse HEAD)" = "$SU2_COMMIT"; \
    cd /src/SU2; \
    ./meson.py setup build -Dwith-mpi=enabled -Dbuildtype=release --prefix=/opt/su2; \
    ./ninja -C build -j"$BUILD_JOBS" install; \
    strip /opt/su2/bin/SU2_CFD; \
    test -x /opt/su2/bin/SU2_CFD; \
    printf '%s %s\n' "$SU2_TAG" "$SU2_COMMIT" > /opt/su2/VERSION

# ---------- стадия 2: образ для расчётов ----------
FROM ${BASE_IMAGE} AS runtime
ARG GIT_SHA=unknown
ENV DEBIAN_FRONTEND=noninteractive
# libgl1, libglu1-mesa, libxrender1, libxcursor1, libxfixes3, libxft2, libfontconfig1, libxinerama1 — нужны
# колесу gmsh (без libgl1 импорт падает: «libGL.so.1: cannot open shared object file», пробная сборка 14.11.2026).
# hadolint ignore=DL3008
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      python3 python3-venv ca-certificates openmpi-bin libopenmpi3t64 calculix-ccx \
      libgl1 libglu1-mesa libxrender1 libxcursor1 libxfixes3 libxft2 libfontconfig1 libxinerama1 libgomp1 \
 && rm -rf /var/lib/apt/lists/*
COPY routes/gtu/platform/docker/requirements.txt /tmp/requirements.txt
RUN --mount=type=secret,id=extra_ca,required=false \
    set -eu; \
    if [ -s /run/secrets/extra_ca ]; then export PIP_CERT=/run/secrets/extra_ca; fi; \
    python3 -m venv /opt/venv; \
    /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt; \
    /opt/venv/bin/python -c "import cantera, numpy, scipy, gmsh, boto3; print('ok', cantera.__version__, gmsh.__version__)"
COPY --from=su2-build /opt/su2 /opt/su2
# Код — пользователю 1000: скрипты трека пишут рядом с собой (run_case.py строит сетку в ../geom/out/).
COPY --chown=1000:1000 routes/gtu /opt/gtu/routes/gtu
# Непривилегированный пользователь: OpenMPI не запускается от root без особых переменных.
# ccx -v печатает версию, но выходит с кодом 201 (проверено 14.11.2026) — проверка без пайпа и без кода выхода.
RUN useradd --create-home --uid 1000 gtu \
 && mkdir -p /work && chown 1000:1000 /work \
 && v="$(ccx -v || true)" && case "$v" in *"Version 2.21"*) ;; *) echo "ccx: $v" >&2; exit 1 ;; esac \
 && printf '%s\n' "$GIT_SHA" > /opt/gtu/GIT_SHA \
 && dpkg-query -W -f '${Package} ${Version}\n' > /opt/gtu/PACKAGES.txt \
 && /opt/venv/bin/pip freeze > /opt/gtu/PIP_FREEZE.txt
LABEL org.opencontainers.image.title="gtu-platform" \
      org.opencontainers.image.revision="${GIT_SHA}" \
      org.opencontainers.image.source="https://github.com/skharkov1246/ocm-research-diary"
# OMPI_MCA_btl_vader_single_copy_mechanism=none: в контейнере без CAP_SYS_PTRACE прямое копирование между
# процессами (CMA) запрещено, обмен идёт через общую память /dev/shm (её размер — sharedMemorySize задания).
ENV PATH=/opt/venv/bin:/opt/su2/bin:$PATH \
    SU2_RUN=/opt/su2/bin \
    GTU_HOME=/opt/gtu \
    GTU_GIT_SHA=${GIT_SHA} \
    OMPI_MCA_btl_vader_single_copy_mechanism=none \
    PYTHONDONTWRITEBYTECODE=1
USER 1000:1000
WORKDIR /work
ENTRYPOINT ["python3", "-u", "/opt/gtu/routes/gtu/platform/batch/entrypoint.py"]
CMD ["--help"]
