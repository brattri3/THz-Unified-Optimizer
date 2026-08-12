# -*- coding: utf-8 -*-
"""tds_track_viewer — интерактивный контроль углового прогона THz-TDS.

Инструмент оператора у спектрометра: показывает угловую зависимость пропускания
ПРЯМО ВО ВРЕМЯ съёмки и раскладывает «дорогу» пустых файлов под следующий прогон.

Это контроль съёмки, а НЕ анализ. Модель не подгоняется, `unified_optimizer` не
импортируется вовсе (числа с ним сверяются, но не берутся у него) — см. ТЗ §9.

ТЗ: `docs/artifacts/C_TZ_tds_track_viewer.md` (зона C, согласовано 2026-08-12).

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m track_viewer.gui        графический интерфейс
    .venv\\Scripts\\python.exe -m track_viewer.cli <путь> печать трека, экспорт CSV
    .venv\\Scripts\\python.exe -m track_viewer.selftest   приёмка §7 ТЗ

Целевая среда — Windows 7 / Python 3.8 / numpy 1.21, поэтому в пакете действуют
ограничения: только `numpy` + `matplotlib` + стандартная библиотека, синтаксис 3.8
(проверяется `tools/check_py38.py`), никаких функций numpy, разошедшихся между
1.x и 2.x (см. `core/compat.py`).
"""
from __future__ import annotations

__version__ = "0.1.0"
