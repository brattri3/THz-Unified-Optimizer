# -*- coding: utf-8 -*-
"""Совместимость numpy 1.21 (целевая машина) ↔ numpy 2.x (машина разработки).

Зачем этот модуль существует
----------------------------
Инструмент разрабатывается на Python 3.14 / numpy 2.5, а работать обязан на
Windows 7 с Python 3.8 / numpy 1.21.6 — последней парой, поддерживаемой этой ОС.
Между этими версиями numpy разошёлся ровно там, где нам нужно:

    np.trapz      есть в 1.21, УДАЛЁН в numpy 2.0
    np.trapezoid  есть в numpy 2.0+, ОТСУТСТВУЕТ в 1.21

То есть любое прямое обращение к интегрированию трапециями ломает инструмент на
одной из двух машин. Ветвление по версии (`getattr(np, 'trapezoid', np.trapz)`)
работает, но прячет проблему в рантайм и даёт разные пути кода на машине, где мы
проверяем числа, и на машине, где владелец их читает. Поэтому здесь — собственная
реализация в три строки: один и тот же код и одинаковые числа везде.
"""
from __future__ import annotations

import io
import sys

import numpy as np


def trapezoid(y, x):
    """Интеграл ``y`` по ``x`` методом трапеций. Замена np.trapz / np.trapezoid.

    Совпадает с обеими версиями numpy бит-в-бит на равномерной сетке: та же
    формула ``Σ (y[i] + y[i+1])/2 · (x[i+1] − x[i])``, то же накопление в float64.
    """
    y = np.asarray(y, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    return float(np.sum((y[1:] + y[:-1]) * 0.5 * np.diff(x)))


def setup_console():
    """Готовит `sys.stdout` к русскому тексту, не ломая консоль Windows.

    Наивное решение «завернуть stdout в UTF-8» **неверно и было причиной
    дефекта**: консоль русской Windows работает в cp866, и UTF-8-байты в ней
    отображаются мусором. Именно поэтому оператору приходилось копировать
    вывод в Блокнот — тот распознаёт UTF-8 и показывает текст правильно.

    Почему обёртку вообще поставили: в выводе есть символы вне cp866 (`Δ`,
    `⚠`, `…`, стрелки), и Python 3.8 на них падает с `UnicodeEncodeError`.
    Правильное лекарство — не менять кодировку, а разрешить **замену**
    неотображаемых символов: русский текст остаётся читаемым, экзотика
    становится `?`.

    Перенаправленный поток (файл, канал) — другое дело: там консоли нет,
    и UTF-8 уместен.
    """
    stream = sys.stdout
    enc = (getattr(stream, "encoding", None) or "").lower()
    if enc in ("utf-8", "utf8"):
        return
    buf = getattr(stream, "buffer", None)
    if buf is None:
        return
    try:
        console = stream.isatty()
    except Exception:                                           # noqa: BLE001
        console = False
    try:
        if console:
            sys.stdout = io.TextIOWrapper(buf, encoding=stream.encoding,
                                          errors="replace", line_buffering=True)
        else:
            sys.stdout = io.TextIOWrapper(buf, encoding="utf-8")
    except Exception:                                           # noqa: BLE001
        pass


class Tee(object):
    """Пишет одновременно в поток и в файл. Нужен для `--log`.

    Лог всегда UTF-8: его открывают Блокнотом, а не консолью, и туда попадают
    все символы без замены. На экране при этом — кодировка консоли.
    """

    def __init__(self, stream, fh):
        self._s = stream
        self._f = fh

    def write(self, text):
        self._s.write(text)
        self._f.write(text)

    def flush(self):
        self._s.flush()
        self._f.flush()


def rfft_freqs(n, dt):
    """Частоты действительного БПФ, ТГц (при ``dt`` в пс).

    Обёртка ради единственного места, где задаётся связь «пс → ТГц»: 1/пс = ТГц,
    поэтому дополнительного множителя нет, и это стоит зафиксировать явно, чтобы
    единицы не пришлось выводить заново при чтении кода.
    """
    return np.fft.rfftfreq(int(n), d=float(dt))
