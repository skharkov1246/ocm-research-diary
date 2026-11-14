# Образ платформы ГТУ — `routes/gtu/platform/docker/`

Один образ для AWS Batch и локального запуска: Ubuntu 24.04, Python 3.12 (venv `/opt/venv`: cantera 3.2.0,
numpy 2.5.3, scipy 1.18.1, gmsh 4.15.2, boto3 1.44.3), CalculiX `ccx` 2.21 (apt), OpenMPI 4.1.6 (apt),
SU2 v8.6.0 «Harrier» с MPI (`/opt/su2`, собирается в первой стадии из тега, коммит 7048fd35 сверяется),
код трека `routes/gtu` → `/opt/gtu/routes/gtu`. Пользователь — uid 1000, рабочая папка `/work`.
Точка входа — `routes/gtu/platform/batch/entrypoint.py` (задание по JSON).

| Файл | Что это |
|---|---|
| `Dockerfile` | две стадии: `su2-build` (сборка SU2) и `runtime` |
| `Dockerfile.dockerignore` | контекст — корень репозитория, в образ идёт только `routes/gtu` без прогонов, сеток и кэшей |
| `requirements.txt` | закреплённые версии пакетов Python |

Почему SU2 собирается, почему Ubuntu 24.04, как собирать за прокси и что проверено — README узла
`../batch/README.md` (раздел «Образ») и `../batch/LOCAL_FIRST.md` (шаги 1–2).

```bash
# из корня чистой копии коммита
docker build -f routes/gtu/platform/docker/Dockerfile \
    --build-arg GIT_SHA="$(git rev-parse HEAD)" -t gtu-platform:"$(git rev-parse --short=12 HEAD)" .
```

В образе для повторяемости: `/opt/gtu/GIT_SHA` (хеш коммита сборки, он же метка `org.opencontainers.image.revision`
и `$GTU_GIT_SHA`), `/opt/su2/VERSION` (тег и коммит SU2), `/opt/gtu/PACKAGES.txt` (версии пакетов apt),
`/opt/gtu/PIP_FREEZE.txt` (версии пакетов Python).
