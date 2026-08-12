# -*- coding: utf-8 -*-
"""Разбор каталога углового прогона: файлы на диске → модель `Scan`.

Правило отбора (ТЗ §2.1)
------------------------
Трассой считается только файл, чьё имя совпало с

    ^(\\d{3})_(a(-?\\d+)|bg)\\.txt$

**Всё, что не совпало, — не трасса, и пропускается молча.** Это не «фильтр
мусора», а сам механизм отсева мета-информации: одним правилом отсекаются
`README.md`, `MEASUREMENT_ORDER.md`, `000_META_TEMPLATE.txt`, `META.txt`, старые
`NNN_META.txt`, подкаталог `microscopy/` — и, что важнее всего, **нормализованные
копии** вида `<набор>_<угол>deg_repN_sig.txt`. Последних в `test_grid_33_11` 108
штук, и они побайтно равны журнальным файлам: попади они в выборку, каждая точка
угловой зависимости удвоилась бы.

Три состояния файла (ТЗ §2.3)
-----------------------------
* **нулевой размер** — «ещё не измерено». Законное состояние: дорога
  раскладывается заранее и заполняется по ходу съёмки. Пропускается молча;
* **разобран** — точка на графике;
* **непустой, но не разобрался** — ошибка. Показывается в панели данных, точки
  на графике нет, прогон не прерывается.

Повторы
-------
`rep` выводится из порядка: первое появление угла среди сигнальных трасс (по
возрастанию сквозного номера) — `rep1`, второе — `rep2`. Без этого повторный
обход тени из блока 4 плана измерений сливается с первым в одну точку, а разница
между ними — это и есть прямая мера дрейфа установки.
"""
from __future__ import annotations

import os
import re
from collections import OrderedDict

import numpy as np

# Сквозной номер + вид. Знак угла — часть имени (`047_a-94.txt`).
TRACE_RE = re.compile(r"^(\d{3})_(?:a(-?\d+)|bg)\.txt$")

EMPTY = "empty"
OK = "ok"
ERROR = "error"


class Trace(object):
    """Одна запись прогона: сигнальная трасса или референс.

    Данные (`t`, `E`) читаются лениво и кэшируются в `Scan`, поэтому объект
    остаётся дешёвым: у прогона их под две сотни, а «Обновить» жмут часто.
    """

    __slots__ = ("seq", "name", "path", "kind", "angle", "rep",
                 "size", "mtime", "state", "error", "t", "E")

    def __init__(self, seq, name, path, kind, angle, size, mtime):
        self.seq = seq
        self.name = name
        self.path = path
        self.kind = kind            # 'sig' | 'bg'
        self.angle = angle          # int для 'sig', None для 'bg'
        self.rep = 1                # проставляется в Scan после сортировки
        self.size = size
        self.mtime = mtime
        self.state = EMPTY if size == 0 else OK   # OK — предварительно, до чтения
        self.error = None
        self.t = None
        self.E = None

    @property
    def key(self):
        """Ключ кэша: имя + размер + время правки (ТЗ §4.4)."""
        return (self.name, self.size, self.mtime)

    @property
    def label(self):
        if self.kind == "bg":
            return u"%s (референс)" % self.name
        if self.rep > 1:
            return u"%s  %+d°  rep%d" % (self.name, self.angle, self.rep)
        return u"%s  %+d°" % (self.name, self.angle)

    def __repr__(self):
        return "<Trace %s %s state=%s>" % (self.name, self.kind, self.state)


def read_trace_file(path):
    """Читает файл трассы → (t, E). Бросает ValueError с внятным текстом.

    Формат: текст, две или три колонки, разделитель — табуляция, шапки нет.
    Третья колонка — позиция линии задержки; решением владельца (O-5, 2026-08-12)
    она читается, но в расчётах не участвует.

    Число точек и длина окна НЕ зашиты: у прогонов они разные (651 точка / 67.7 пс
    у `test_grid_33_11`, но окно настраивается на приборе), поэтому берутся из
    файла. `np.loadtxt` сам разбирает и табуляцию, и пробелы.
    """
    try:
        arr = np.loadtxt(path, dtype=np.float64, ndmin=2)
    except Exception as exc:                                    # noqa: BLE001
        # Текст numpy подставляем только если он ASCII: он цитирует содержимое
        # файла, а нечисловой файл здесь обычно русский текст в cp1251, и в
        # панели данных это выглядит как «СЌС‚Рѕ» — хуже, чем без подробностей.
        detail = u"%s" % exc
        try:
            detail.encode("ascii")
        except UnicodeEncodeError:
            detail = u"файл содержит не числа"
        raise ValueError(u"не разобрался как числовая таблица: %s" % detail)

    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(u"ожидались минимум 2 колонки, найдено %s"
                         % (arr.shape[1] if arr.ndim == 2 else u"?"))
    if arr.shape[0] < 8:
        raise ValueError(u"слишком короткая трасса: %d точек" % arr.shape[0])

    t = arr[:, 0]
    E = arr[:, 1]
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(E)):
        raise ValueError(u"в файле есть NaN или бесконечность")
    if np.any(np.diff(t) <= 0):
        raise ValueError(u"колонка времени не возрастает монотонно")
    return t, E


class Scan(object):
    """Каталог одного прогона.

    Хранит список трасс и кэш прочитанных данных. `refresh()` перечитывает
    каталог: новые и изменившиеся файлы читаются с диска, прежние берутся из
    кэша по ключу (имя, размер, mtime).

    Автообновления по таймеру здесь нет и не будет без отдельной просьбы
    владельца (ТЗ §4.4): чтение каталога в момент, когда прибор дописывает
    трассу, вернёт обрезанный файл, и точка молча уедет.
    """

    def __init__(self, directory):
        self.directory = str(directory)
        self.traces = []            # все трассы, по возрастанию seq
        # key -> (True, (t, E)) | (False, текст ошибки). Флаг отдельным элементом,
        # а не меткой в первом поле: там лежит массив, и `==` дал бы поэлементное
        # сравнение вместо проверки признака.
        self._cache = OrderedDict()
        self.refresh()

    # ------------------------------------------------------------------ чтение
    def refresh(self, drop_cache=False):
        """Перечитать каталог. `drop_cache=True` — полный пересчёт (ТЗ §4.4)."""
        if drop_cache:
            self._cache.clear()

        try:
            names = os.listdir(self.directory)
        except OSError as exc:
            raise ValueError(u"каталог не читается: %s" % exc)

        found = []
        for name in names:
            m = TRACE_RE.match(name)
            if m is None:
                continue                    # не трасса — молча мимо (ТЗ §2.1)
            path = os.path.join(self.directory, name)
            if not os.path.isfile(path):
                continue
            st = os.stat(path)
            seq = int(m.group(1))
            if m.group(2) is None:
                found.append(Trace(seq, name, path, "bg", None, st.st_size, st.st_mtime))
            else:
                found.append(Trace(seq, name, path, "sig", int(m.group(2)),
                                   st.st_size, st.st_mtime))

        found.sort(key=lambda tr: (tr.seq, tr.name))
        self._assign_reps(found)
        for tr in found:
            self._load(tr)
        self.traces = found
        return self

    def _assign_reps(self, traces):
        """rep = порядковый номер появления угла среди сигнальных трасс."""
        seen = {}
        for tr in traces:
            if tr.kind != "sig":
                continue
            seen[tr.angle] = seen.get(tr.angle, 0) + 1
            tr.rep = seen[tr.angle]

    def _load(self, tr):
        if tr.size == 0:
            tr.state = EMPTY
            return
        cached = self._cache.get(tr.key)
        if cached is None:
            try:
                cached = (True, read_trace_file(tr.path))
            except ValueError as exc:
                cached = (False, u"%s" % exc)
            self._cache[tr.key] = cached
        if cached[0]:
            tr.state = OK
            tr.error = None
            tr.t, tr.E = cached[1]
        else:
            tr.state = ERROR
            tr.error = cached[1]
            tr.t = tr.E = None

    # ------------------------------------------------------------------ выборки
    @property
    def signals(self):
        """Сигнальные трассы в состоянии OK — то, что перебирает слайдер (ТЗ §4.2)."""
        return [tr for tr in self.traces if tr.kind == "sig" and tr.state == OK]

    @property
    def backgrounds(self):
        return [tr for tr in self.traces if tr.kind == "bg" and tr.state == OK]

    @property
    def errors(self):
        return [tr for tr in self.traces if tr.state == ERROR]

    def counts(self):
        """Сводка для панели «N/M снято»."""
        sig = [tr for tr in self.traces if tr.kind == "sig"]
        bg = [tr for tr in self.traces if tr.kind == "bg"]
        return {
            "sig_total": len(sig),
            "sig_done": len([tr for tr in sig if tr.state == OK]),
            "bg_total": len(bg),
            "bg_done": len([tr for tr in bg if tr.state == OK]),
            "empty": len([tr for tr in self.traces if tr.state == EMPTY]),
            "errors": len(self.errors),
        }

    def neighbours_bg(self, seq):
        """Ближайшие референсы по обе стороны от номера: (ранний, поздний).

        Любой из двух может быть None — у первой трассы прогона нет раннего
        референса, у последней нет позднего. Обработка этого случая — на стороне
        `physics.transmission`, и она обязана пометить откат явно (ТЗ §3.5).
        """
        earlier = later = None
        for tr in self.backgrounds:
            if tr.seq < seq:
                if earlier is None or tr.seq > earlier.seq:
                    earlier = tr
            elif tr.seq > seq:
                if later is None or tr.seq < later.seq:
                    later = tr
        return earlier, later

    def next_unmeasured(self):
        """Первый по порядку файл нулевого размера — «следующее действие» оператора."""
        for tr in self.traces:
            if tr.state == EMPTY:
                return tr
        return None
