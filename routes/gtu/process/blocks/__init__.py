"""Библиотека блоков процесса «от топлива до вала» двигателя-заменителя (трек «Простая литая ГТУ-заменитель», issue #59).

Блоки: fuel (подготовка топлива), compressor, combustor, turbine, hx (теплообменники), shaft (вал и стык пакета);
совместимость - compat; цепочка - chain; единый вызов - evaluate(concept) (evaluator.py). Уровень облика, всё - оценка.
"""
import os as _os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):     # одна нить BLAS: на общей машине
    _os.environ.setdefault(_v, "1")                                         # cycles.lam() иначе идёт 16 с вместо 0,05

from .concept import DEFAULT, normalize, label  # noqa: F401
from .evaluator import evaluate  # noqa: F401
