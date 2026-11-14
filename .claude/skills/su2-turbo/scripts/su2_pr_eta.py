"""Полная степень сжатия, изоэнтропический КПД и дисбаланс расхода по history.csv SU2.

Запуск: python3 su2_pr_eta.py history.csv [gamma] [маркер_входа] [маркер_выхода]
Нужны в HISTORY_OUTPUT группы FLOW_COEFF_SURF (или поля SURFACE_MASSFLOW, SURFACE_TOTAL_PRESSURE,
SURFACE_TOTAL_TEMPERATURE) и MARKER_ANALYZE= (вход, выход). Берётся последняя строка файла.
Формула КПД — для идеального газа с постоянным gamma (как FLUID_MODEL= IDEAL_GAS).
Нормаль границы в SU2 направлена внутрь области: при прямом течении расход на входе > 0, на выходе < 0.
Проверено 14.11.2026 на history.csv SU2 v8.6.0 (пробный прогон колеса core130, coarse).
"""
import sys

path = sys.argv[1]
g = float(sys.argv[2]) if len(sys.argv) > 2 else 1.4
mi = sys.argv[3] if len(sys.argv) > 3 else "inlet"
mo = sys.argv[4] if len(sys.argv) > 4 else "outlet"

lines = [ln for ln in open(path, encoding="utf-8") if ln.strip()]
head = [h.strip().strip('"') for h in lines[0].split(",")]
row = dict(zip(head, (float(v) for v in lines[-1].split(","))))

pt1, pt2 = row[f"Avg_TotalPress({mi})"], row[f"Avg_TotalPress({mo})"]
tt1, tt2 = row[f"Avg_TotalTemp({mi})"], row[f"Avg_TotalTemp({mo})"]
m1, m2 = row[f"Avg_Massflow({mi})"], row[f"Avg_Massflow({mo})"]

pr = pt2 / pt1
eta = (pr ** ((g - 1) / g) - 1) / (tt2 / tt1 - 1)
print(f"PR_tt {pr:.4f}  eta_tt {eta:.4f}  m_in {m1:+.4f}  m_out {m2:+.4f}  "
      f"дисбаланс {abs(m1 + m2) / abs(m1):.3%}")
if not (m1 > 0 > m2):
    print("ВНИМАНИЕ: знаки расхода не «вход > 0, выход < 0» — обратное течение или расчёт не установился; "
          "PR и eta не публиковать")
