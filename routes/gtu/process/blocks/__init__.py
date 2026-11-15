"""Библиотека блоков процесса «от топлива до вала» двигателя-заменителя (трек «Простая литая ГТУ-заменитель», issue #59).

Блоки: fuel (подготовка топлива), compressor, combustor, turbine, hx (теплообменники), shaft (вал и стык пакета);
совместимость - compat; цепочка - chain; единый вызов - evaluate.evaluate(concept). Уровень облика, всё - оценка.
"""
from .concept import DEFAULT, normalize, label  # noqa: F401
from .evaluate import evaluate  # noqa: F401
