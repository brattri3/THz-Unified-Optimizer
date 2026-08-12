# -*- coding: utf-8 -*-
"""Сборка офлайн-комплекта track_viewer для Windows 7.

Почему сборка скриптом, а не руками
-----------------------------------
`win7_attenuator_kit` собирали вручную, и он разошёлся с репозиторием: остался
на старых именах наборов, а `.exe` в нём собран под Python 3.14 и на Windows 7
не запускается в принципе. Ручной комплект нельзя пересобрать после правки кода,
поэтому он устаревает молча. Этот скрипт можно запустить заново в любой момент.

Что делает
----------
1. копирует пакет `track_viewer` (без `__pycache__` и без служебных инструментов);
2. копирует данные двух эталонных прогонов БЕЗ микрофото — чтобы главный пункт
   приёмки (T = 98.9089 %) проверялся на самой машине у прибора, а не только
   на машине разработчика;
3. кладёт колёса под `cp38` для обеих разрядностей (скачивает или берёт из кэша);
4. переиспользует установщики Python 3.8.10 из `win7_attenuator_kit`, если тот
   рядом: это те же самые файлы, качать их второй раз незачем;
5. пишет `.bat`-запускалки и `README_WIN7.txt`.

Требуется интернет — **на этой машине и заранее**, а не у прибора. Если колёс
под 3.8 не окажется, узнать об этом надо здесь.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe track_viewer/tools/make_win7_kit.py
    .venv\\Scripts\\python.exe track_viewer/tools/make_win7_kit.py --out <путь>
    .venv\\Scripts\\python.exe track_viewer/tools/make_win7_kit.py --no-download
"""
from __future__ import annotations

import argparse
import io
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))    # track_viewer/
REPO = os.path.dirname(HERE)
PARENT = os.path.dirname(REPO)

NUMPY = "numpy==1.21.6"          # последняя, поддерживающая Python 3.8
MATPLOTLIB = "matplotlib==3.5.3"  # последняя, поддерживающая Python 3.8
PLATFORMS = ("win_amd64", "win32")

DATASETS = ("test_grid_33_11", "purewave")

# Файлы пакета, которые в поставку не идут: это инструменты разработчика.
SKIP_FILES = ("smoke_gui.py", "make_win7_kit.py")

README = u"""\
=============================================================================
  ОФЛАЙН-КОМПЛЕКТ  tds_track_viewer  ДЛЯ WINDOWS 7
  контроль углового прогона THz-TDS прямо во время съёмки
=============================================================================

  Интернет не нужен: в комплекте есть всё, включая установщик Python.

  Почему не .exe. Единственный Python, работающий на Windows 7, — 3.8;
  собранный на более новом .exe там не запустится в принципе (именно это
  случилось с комплектом аттенюатора v0.3). Отдельный плюс комплекта: код
  правится блокнотом без пересборки, и антивирус к нему не придирается.

-----------------------------------------------------------------------------
 СОСТАВ
-----------------------------------------------------------------------------
  python-3.8.10-amd64.exe        установщик Python, 64 бита
  python-3.8.10-win32.exe        установщик Python, 32 бита
  wheels\\win_amd64\\               библиотеки под 64 бита (10 файлов)
  wheels\\win32\\                   библиотеки под 32 бита (10 файлов)
  app\\track_viewer\\               сама программа
  app\\data_pool\\                  два эталонных прогона для самопроверки
  install.bat                    установка библиотек (один раз)
  start.bat                      запуск программы
  selftest.bat                   самопроверка
  console.bat                    консольный режим и выгрузка CSV

-----------------------------------------------------------------------------
 УСТАНОВКА (один раз, примерно 5 минут)
-----------------------------------------------------------------------------
 1. Узнать разрядность Windows: Пуск -> Компьютер -> правой кнопкой ->
    Свойства -> «Тип системы». Нужна Windows 7 SP1.

 2. Запустить ПОДХОДЯЩИЙ установщик:
      64 бита -> python-3.8.10-amd64.exe
      32 бита -> python-3.8.10-win32.exe
    ВАЖНО: внизу поставить галочку [x] Add Python 3.8 to PATH,
    затем «Install Now».

 3. Двойной щелчок по install.bat — он сам определит разрядность и поставит
    библиотеки из папки wheels. Интернет не потребуется.

-----------------------------------------------------------------------------
 ПРОВЕРКА ПОСЛЕ УСТАНОВКИ
-----------------------------------------------------------------------------
  Двойной щелчок по selftest.bat. Ожидается строка

      === пройдено 27 из 27 ===

  Главное, что там проверяется: пропускание в ярком положении набора
  test_grid_33_11 равно 98.9089 %. Ровно это же число получено на машине
  разработчика и сходится с расчётом в ядре проекта. Если оно другое —
  дальше не идти, разбираться.

-----------------------------------------------------------------------------
 РАБОТА
-----------------------------------------------------------------------------
  start.bat        открывает окно программы.
                   Кнопка «выбрать…» — указать каталог текущего прогона.
                   Кнопка «ОБНОВИТЬ» — перечитать каталог после новых трасс.
                   Кнопка «разложить дорогу…» — создать пустые файлы под
                   следующий прогон.

  ОБНОВЛЕНИЕ ПО ТАЙМЕРУ НЕ ДЕЛАЕТСЯ НАМЕРЕННО: чтение каталога в момент,
  когда прибор дописывает трассу, вернуло бы обрезанный файл. Жать
  «ОБНОВИТЬ» вручную, когда трасса дописана.

  console.bat      консольный режим: печать всего трека и выгрузка в CSV
                   для Origin Pro. Пример:
                     python -m track_viewer.cli C:\\путь\\к\\прогону --csv трек.csv

-----------------------------------------------------------------------------
 ЧТО ПРОГРАММА НЕ ДЕЛАЕТ
-----------------------------------------------------------------------------
  * не подгоняет модель — это контроль съёмки, а не анализ;
  * не изменяет и не удаляет измеренные трассы ни при каких настройках;
  * пишет в каталог данных только пустые файлы, META.txt и README.md, и
    только при явной генерации дороги;
  * не заменяет предполётный чек-лист: центровку образца проверяют руками,
    программа этого не видит.

-----------------------------------------------------------------------------
 ЕСЛИ ЧТО-ТО НЕ ТАК
-----------------------------------------------------------------------------
  «python не является внутренней или внешней командой» — при установке не
  поставили галочку «Add Python 3.8 to PATH». Переустановить с галочкой.

  install.bat ругается на разрядность — проверить, что установлен Python
  той же разрядности, что и папка wheels, из которой он ставит.

  Окно открывается, но графиков нет — запустить selftest.bat и прислать
  его вывод целиком.
"""

INSTALL_BAT = """\
@echo off
chcp 1251 >nul
cd /d "%~dp0"
echo.
echo Определяю разрядность Python...
python -c "import struct,sys; sys.stdout.write('win_amd64' if struct.calcsize('P')==8 else 'win32')" > "%TEMP%\\tv_arch.txt" 2>nul
if errorlevel 1 goto nopython
set /p ARCH=<"%TEMP%\\tv_arch.txt"
del "%TEMP%\\tv_arch.txt"
echo Python: %ARCH%
echo.
echo Ставлю библиотеки из wheels\\%ARCH% ...
python -m pip install --no-index --find-links "wheels\\%ARCH%" numpy matplotlib
if errorlevel 1 goto failed
echo.
echo Готово. Теперь запустите selftest.bat
pause
exit /b 0

:nopython
echo.
echo ОШИБКА: python не найден.
echo При установке Python нужно было поставить галочку "Add Python 3.8 to PATH".
echo Переустановите Python 3.8.10 из этой папки с этой галочкой.
pause
exit /b 1

:failed
echo.
echo ОШИБКА при установке библиотек. Вывод выше.
pause
exit /b 1
"""

START_BAT = """\
@echo off
chcp 1251 >nul
cd /d "%~dp0app"
python -m track_viewer.gui
if errorlevel 1 pause
"""

SELFTEST_BAT = """\
@echo off
chcp 1251 >nul
cd /d "%~dp0app"
set PYTHONIOENCODING=utf-8
python -m track_viewer.selftest
echo.
pause
"""

CONSOLE_BAT = """\
@echo off
chcp 1251 >nul
cd /d "%~dp0app"
set PYTHONIOENCODING=utf-8
echo Консольный режим track_viewer.
echo.
echo   python -m track_viewer.cli data_pool\\test_grid_33_11
echo   python -m track_viewer.cli ПУТЬ --csv трек.csv
echo   python -m track_viewer.cli ПУТЬ --dc none --ref average
echo.
cmd /k
"""


def copy_package(dst):
    """Копирует пакет без кэша и без инструментов разработчика."""
    n = 0
    for dirpath, dirnames, filenames in os.walk(HERE):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        rel = os.path.relpath(dirpath, HERE)
        out = os.path.join(dst, "track_viewer") if rel == "." else \
            os.path.join(dst, "track_viewer", rel)
        if not os.path.isdir(out):
            os.makedirs(out)
        for name in filenames:
            if not name.endswith(".py") or name in SKIP_FILES:
                continue
            shutil.copy2(os.path.join(dirpath, name), os.path.join(out, name))
            n += 1
    return n


def copy_data(dst):
    """Эталонные прогоны без микрофото: 6 МБ против 125 МБ с ними."""
    n = 0
    for ds in DATASETS:
        src = os.path.join(REPO, "data_pool", ds)
        if not os.path.isdir(src):
            print(u"  ПРЕДУПРЕЖДЕНИЕ: нет %s — самопроверка на месте не пройдёт" % src)
            continue
        out = os.path.join(dst, "data_pool", ds)
        if not os.path.isdir(out):
            os.makedirs(out)
        for name in sorted(os.listdir(src)):
            path = os.path.join(src, name)
            if not os.path.isfile(path):
                continue            # microscopy/ и прочие подкаталоги — мимо
            shutil.copy2(path, os.path.join(out, name))
            n += 1
    return n


def fetch_wheels(dst, download=True):
    total = 0
    for platform in PLATFORMS:
        out = os.path.join(dst, platform)
        if not os.path.isdir(out):
            os.makedirs(out)
        have = [n for n in os.listdir(out) if n.endswith(".whl")]
        if have and not download:
            print(u"  %s: уже лежит %d колёс, скачивание пропущено"
                  % (platform, len(have)))
            total += len(have)
            continue
        cmd = [sys.executable, "-m", "pip", "download",
               "--only-binary=:all:", "--python-version", "3.8",
               "--implementation", "cp", "--abi", "cp38",
               "--platform", platform, "-d", out, NUMPY, MATPLOTLIB]
        print(u"  %s: качаю…" % platform)
        res = subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if res != 0:
            raise RuntimeError(
                u"не удалось скачать колёса под %s. Это надо решить ЗДЕСЬ, "
                u"на машине с интернетом, а не у прибора." % platform)
        got = [n for n in os.listdir(out) if n.endswith(".whl")]
        print(u"  %s: %d колёс" % (platform, len(got)))
        total += len(got)
    return total


def copy_python_installers(dst):
    """Берёт установщики из комплекта аттенюатора, если он рядом."""
    src = os.path.join(PARENT, "win7_attenuator_kit")
    found = []
    if os.path.isdir(src):
        for name in os.listdir(src):
            if re.match(r"^python-3\.8\.\d+-(amd64|win32)\.exe$", name):
                target = os.path.join(dst, name)
                if not os.path.exists(target):
                    shutil.copy2(os.path.join(src, name), target)
                found.append(name)
    return found


def write_text(path, text):
    with io.open(path, "w", encoding="cp1251", errors="replace",
                 newline="\r\n") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(PARENT, "win7_track_viewer_kit"))
    ap.add_argument("--no-download", action="store_true",
                    help=u"не качать колёса, если они уже лежат в комплекте")
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    print(u"комплект: %s" % out)
    if not os.path.isdir(out):
        os.makedirs(out)

    app = os.path.join(out, "app")
    print(u"пакет:")
    print(u"  скопировано файлов: %d" % copy_package(app))
    print(u"данные для самопроверки:")
    print(u"  скопировано файлов: %d (микрофото исключены)" % copy_data(app))
    print(u"колёса:")
    print(u"  всего: %d" % fetch_wheels(os.path.join(out, "wheels"),
                                        download=not args.no_download))
    print(u"установщики Python:")
    inst = copy_python_installers(out)
    if inst:
        for name in inst:
            print(u"  %s" % name)
    else:
        print(u"  НЕ НАЙДЕНЫ. Положите python-3.8.10-amd64.exe и "
              u"python-3.8.10-win32.exe в корень комплекта вручную:")
        print(u"  https://www.python.org/downloads/release/python-3810/")

    # Тексты для Windows 7: cp1251 и CRLF — «Блокнот» там не понимает UTF-8
    # без BOM и склеивает строки с одними LF.
    write_text(os.path.join(out, "README_WIN7.txt"), README)
    write_text(os.path.join(out, "install.bat"), INSTALL_BAT)
    write_text(os.path.join(out, "start.bat"), START_BAT)
    write_text(os.path.join(out, "selftest.bat"), SELFTEST_BAT)
    write_text(os.path.join(out, "console.bat"), CONSOLE_BAT)
    print(u"запускалки и README_WIN7.txt: записаны (cp1251, CRLF)")

    # Кэш байт-кода мог остаться от прогона приёмки внутри комплекта: он собран
    # текущим Python, а на целевой машине 3.8 — там он не подойдёт и только
    # запутает. Пересборка комплекта обязана его убирать.
    dropped = 0
    for dirpath, dirnames, _filenames in os.walk(out):
        for name in list(dirnames):
            if name == "__pycache__":
                shutil.rmtree(os.path.join(dirpath, name), ignore_errors=True)
                dirnames.remove(name)
                dropped += 1
    if dropped:
        print(u"кэш байт-кода: удалено каталогов %d" % dropped)

    size = 0
    for dirpath, _dirnames, filenames in os.walk(out):
        for name in filenames:
            size += os.path.getsize(os.path.join(dirpath, name))
    print(u"")
    print(u"готово, объём комплекта: %.0f МБ" % (size / 1024.0 / 1024.0))
    print(u"перенести целиком на ПК спектрометра и открыть README_WIN7.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
