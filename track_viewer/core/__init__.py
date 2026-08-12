# -*- coding: utf-8 -*-
"""Ядро track_viewer: разбор каталога прогона и расчёты.

**Инвариант пакета: здесь нет ни `import tkinter`, ни `import matplotlib`.**
Ровно это позволяет прогнать всю приёмку из `selftest.py` без открытия окна и
без графического окружения вообще. Нарушение инварианта проверяется тестом
`test_core_is_headless` в `selftest.py` — не «по договорённости», а машинно.
"""
from __future__ import annotations
