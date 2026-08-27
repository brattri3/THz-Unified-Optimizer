# -*- coding: utf-8 -*-
"""Палитра и оформление окна.

Цвета берутся из `attenuator_app.core.plots` -- того же модуля, которым
рисуются PNG в остальном проекте. Один источник на два рендерера: иначе экран
и печатная картинка разъедутся молча при первой же правке палитры.

Импорт `core.plots` НЕ тянет matplotlib: там он подгружается лениво, а сами
цвета лежат обычными строками на уровне модуля.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.core import plots                     # noqa: E402

SERIES = plots.SERIES          # синий / оранжевый / бирюзовый
SURFACE = plots.SURFACE        # фон поля графика
INK = plots.INK                # основной текст
INK2 = plots.INK2              # вторичный текст
MUTED = plots.MUTED            # подписи осей
GRID = plots.GRID              # сетка
AXIS = plots.AXIS              # рамки и линии осей
STATUS = plots.STATUS          # good / warning / critical

MODEL = SERIES[0]              # кривая прибора
MARK = SERIES[1]               # точка запроса
PANEL = "#f4f3ef"              # фон панелей окна
LINE = "#dcdbd4"               # разделители

#: запас над 100 % и над 0 дБ. При наклонённом источнике максимум уезжает с
#: нуля шкал, и кривая слегка превышает опорный отсчёт -- это верно физически,
#: и обрезать макушку нельзя (решение владельца 2026-08-27)
DB_LIMITS = (-55.0, 1.0)
PCT_LIMITS = (0.0, 105.0)

FONT_UI = "Segoe UI"
FONT_MONO = "Cascadia Mono"

QSS = """
QWidget { background: %(panel)s; color: %(ink)s; font-family: '%(font)s'; font-size: 12px; }
QGroupBox { border: 1px solid %(line)s; border-radius: 3px; margin-top: 14px;
            padding: 6px 6px 8px 6px; background: %(surface)s; }
QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 3px;
                   color: %(ink2)s; font-size: 11px; font-weight: 600; }
QLabel#hint { color: %(muted)s; font-size: 11px; }
QLabel#warn { color: %(warn)s; font-size: 11px; }
QDoubleSpinBox, QLineEdit { background: #ffffff; border: 1px solid %(axis)s;
                            border-radius: 2px; padding: 3px 4px;
                            font-family: '%(mono)s'; }
QDoubleSpinBox:disabled, QLabel:disabled, QRadioButton:disabled { color: %(muted)s; }
QPushButton { background: %(panel)s; border: 1px solid %(axis)s; border-radius: 2px;
              padding: 5px 12px; }
QPushButton:hover { background: #ecebe6; }
QScrollArea { border: none; background: %(panel)s; }
QStatusBar { background: %(panel)s; border-top: 1px solid %(line)s; color: %(ink2)s; }
""" % {"panel": PANEL, "ink": INK, "ink2": INK2, "muted": MUTED, "line": LINE,
       "surface": SURFACE, "axis": AXIS, "font": FONT_UI, "mono": FONT_MONO,
       "warn": STATUS["warning"]}


def configure_pyqtgraph() -> None:
    """Глобальные настройки pyqtgraph под палитру проекта."""
    import pyqtgraph as pg

    pg.setConfigOption("background", SURFACE)
    pg.setConfigOption("foreground", INK2)
    pg.setConfigOption("antialias", True)


def pens():
    """Перья графика. Отдельной функцией -- Qt должен быть уже поднят."""
    import pyqtgraph as pg
    from PySide6.QtCore import Qt

    return {
        "model": pg.mkPen(MODEL, width=2),
        "grid": pg.mkPen(GRID, width=1),
        "cursor": pg.mkPen(INK2, width=1, style=Qt.DashLine),
        "mark": pg.mkBrush(MARK),
    }


def limits(units: str) -> tuple[float, float]:
    return DB_LIMITS if units == "dB" else PCT_LIMITS


def axis_label(units: str) -> str:
    return "attenuation, dB" if units == "dB" else "transmission, %"
