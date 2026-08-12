# -*- coding: utf-8 -*-
"""Проверка синтаксической совместимости пакета с Python 3.8.

Зачем
-----
Инструмент разрабатывается на Python 3.14, а работать обязан на Windows 7, где
последний поддерживаемый Python — 3.8. Синтаксис, появившийся после 3.8
(`match`, `x: int = 1` в позиционно-только форме, `except*`, оператор `|` для
типов в рантайме и прочее), пройдёт на машине разработчика и упадёт у прибора —
то есть ровно там, где отладка дороже всего.

`ast.parse(..., feature_version=(3, 8))` заставляет парсер CPython отвергать
конструкции новее 3.8, не запуская код. Это ловит синтаксис; **семантику
библиотек он не ловит** (например, отсутствие `np.trapezoid` в numpy 1.21) —
для этого есть `core/compat.py` и приёмка на самой машине.

История: заявление «проверено на совместимость с 3.8, 16 файлов, 0 ошибок» в
`win7_attenuator_kit/README_WIN7.txt` было разовой ручной проверкой; инструмента,
который можно перезапустить после правки, в зоне C до сих пор не существовало.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe track_viewer/tools/check_py38.py
    .venv\\Scripts\\python.exe track_viewer/tools/check_py38.py attenuator_app
"""
from __future__ import annotations

import ast
import io
import os
import sys

TARGET = (3, 8)


def check_file(path):
    """-> None если файл совместим, иначе строка с описанием ошибки."""
    with io.open(path, "r", encoding="utf-8") as fh:
        src = fh.read()
    try:
        ast.parse(src, filename=path, feature_version=TARGET)
    except SyntaxError as exc:
        return u"строка %s: %s" % (exc.lineno, exc.msg)
    return None


def walk(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in ("__pycache__", "build", "dist", ".git")]
        for name in sorted(filenames):
            if name.endswith(".py"):
                yield os.path.join(dirpath, name)


def main(argv):
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    roots = argv[1:] or [here]

    checked = 0
    bad = []
    for root in roots:
        if not os.path.exists(root):
            print(u"нет такого пути: %s" % root)
            return 2
        for path in walk(root):
            checked += 1
            err = check_file(path)
            if err:
                bad.append((path, err))

    print(u"проверено файлов: %d (цель — Python %d.%d)" % (checked, TARGET[0], TARGET[1]))
    if bad:
        for path, err in bad:
            print(u"  НЕСОВМЕСТИМО  %s — %s" % (os.path.relpath(path, here), err))
        print(u"итог: %d несовместимых" % len(bad))
        return 1
    print(u"итог: несовместимого синтаксиса нет")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
