# -*- coding: utf-8 -*-
"""Узкая точка входа для запускалок и будущей сборки.

Существует по той же причине, что и `attenuator_app/gui_entry.py`: `.bat`-файл
комплекта и сборщик исполняемого файла должны указывать на модуль, который
можно запустить и как скрипт (`python gui_entry.py`), и как модуль пакета
(`python -m track_viewer.gui_entry`). Сам `gui.py` для этого не годится —
относительные импорты (`from .core import ...`) при прямом запуске файла падают.
"""
from __future__ import annotations

import os
import sys


def main():
    # Прямой запуск файла: корень пакета в sys.path, дальше обычный импорт.
    if __package__ in (None, ""):
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from track_viewer.gui import main as gui_main
    else:
        from .gui import main as gui_main
    return gui_main()


if __name__ == "__main__":
    sys.exit(main())
