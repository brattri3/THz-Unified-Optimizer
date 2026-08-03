#!/usr/bin/env python
"""Preflight-проверка-и-занятие лока роли (CHARTER §12), через файл-мьютекс.

Версия 2 (2026-08-03, песочница coord-debug): реальный тест показал, что проверка по
markdown-полю "Статус" в sessions/<ID>.md — это read-then-write без атомарности: два
процесса читают FREE до того, как любой из них успевает записать, и оба проходят гейт
(lost update). Заменено на `coordination/sessions/<ID>.lock` — создаётся ТОЛЬКО через
эксклюзивное атомарное открытие (os.O_CREAT | os.O_EXCL), которое ОС гарантированно
проваливает для второго создателя независимо от гонки по времени.

Usage:
    python coordination/tools/preflight.py <ROLE_ID>                  # только проверить
    python coordination/tools/preflight.py <ROLE_ID> --acquire <NICK> # проверить-и-занять

Коды возврата:
    0 — FREE (проверка) / лок только что занят с нуля (занятие)
    1 — лок просрочен: (проверка) реклеймабелен, ничего не менялось /
        (занятие) был просрочен, реклеймлен ЭТИМ вызовом -> ОБЯЗАТЕЛЬНО [STALE-LOCK] в ACTIVITY.md
    2 — ЗАНЯТО: держатель свежий (< таймаута) -> НЕ стартовать. В режиме занятия файл
        .lock НЕ создан и не тронут.
    3 — неопределённое состояние (испорченный .lock, либо гонка ровно в момент реклейма
        удалила и пересоздала лок раньше нас) -> решение за человеком, ничего не трогать
    4 — ошибка использования / файл сессии не найден
"""
import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SESSIONS_DIR = REPO_ROOT / "coordination" / "sessions"
NICKS_FILE = REPO_ROOT / "coordination" / "NICKS.md"

TS_FORMAT = "%Y-%m-%d %H:%M:%S"
DEFAULT_TIMEOUT_MIN = 15


def lock_path(role: str) -> Path:
    return SESSIONS_DIR / f"{role}.lock"


def read_lock(path: Path):
    """Возвращает (nick, datetime) или (None, None) если файл пуст/не разбирается."""
    try:
        text = path.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        return None, None
    lines = text.splitlines()
    if len(lines) < 2:
        return None, None
    nick = lines[0].strip()
    try:
        ts = datetime.strptime(lines[1].strip(), TS_FORMAT)
    except ValueError:
        return nick, None
    return nick, ts


def write_lock_atomic(path: Path, nick: str, now: datetime) -> bool:
    """True если этот вызов реально СОЗДАЛ файл (выиграл гонку); False если файл уже существовал."""
    try:
        fd = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(f"{nick}\n{now.strftime(TS_FORMAT)}\n")
    return True


def check_nicks_registered(nick: str) -> bool:
    if not NICKS_FILE.exists() or not nick:
        return True
    return nick in NICKS_FILE.read_text(encoding="utf-8")


def inspect(role: str, timeout_min: int) -> int:
    path = lock_path(role)
    if not path.exists():
        print(f"[preflight] Роль {role}: .lock отсутствует -> FREE.")
        return 0
    nick, ts = read_lock(path)
    if ts is None:
        print(f"[preflight] Роль {role}: .lock существует, но не разобрать держателя/время "
              f"(содержимое повреждено?) -> решение за человеком.")
        return 3
    age_min = (datetime.now() - ts).total_seconds() / 60.0
    print(f"[preflight] Роль {role}: держатель {nick!r}, взят {ts.strftime(TS_FORMAT)} "
          f"({age_min:.1f} мин назад, таймаут {timeout_min} мин)")
    if age_min < timeout_min:
        print("[preflight] РЕЗУЛЬТАТ: ЗАНЯТО (свежий лок). Не стартовать.")
        return 2
    print("[preflight] РЕЗУЛЬТАТ: ПРОСРОЧЕНО. Реклейм разрешён (но не выполнен — это read-only режим).")
    return 1


def acquire(role: str, nick: str, timeout_min: int) -> int:
    path = lock_path(role)
    now = datetime.now()

    if write_lock_atomic(path, nick, now):
        print(f"[preflight] Роль {role}: .lock не существовал -> СОЗДАН с нуля, держатель {nick}.")
        print("[preflight] РЕЗУЛЬТАТ: ЗАНЯТО ЭТИМ ВЫЗОВОМ (0). Можно приступать к работе под ролью.")
        return 0

    holder, ts = read_lock(path)
    if ts is None:
        print(f"[preflight] Роль {role}: .lock существует, но содержимое не разбирается "
              f"(держатель {holder!r}) -> НЕ трогать, решение за человеком.")
        return 3

    age_min = (now - ts).total_seconds() / 60.0
    print(f"[preflight] Роль {role}: .lock уже занят -> держатель {holder!r}, взят "
          f"{ts.strftime(TS_FORMAT)} ({age_min:.1f} мин назад, таймаут {timeout_min} мин).")

    if age_min < timeout_min:
        print("[preflight] РЕЗУЛЬТАТ: ЗАНЯТО (2). НЕ ЗАНИМАТЬ, не создавать файл, не стартовать работу.")
        return 2

    # Просрочено -> реклейм: удалить и создать заново. Небольшое остаточное окно гонки
    # между remove() и повторным open(..., O_EXCL) намеренно не устраняется транзакционно
    # (санкция владельца, песочница coord-debug 2026-08-03) - при коллизии ровно в этом
    # окне ниже сработает ветка "гонка при реклейме" и вернёт код 3.
    try:
        path.unlink()
    except FileNotFoundError:
        pass  # кто-то другой уже убрал - ниже решит повторный os.open

    if write_lock_atomic(path, nick, now):
        print(f"[preflight] РЕЗУЛЬТАТ: РЕКЛЕЙМ (1). Прежний держатель {holder!r} был просрочен, "
              f"лок пересоздан на {nick}. ОБЯЗАТЕЛЬНО запиши [STALE-LOCK] в ACTIVITY.md + "
              "хэндофф к ORCH в HANDOFFS.md, ПОТОМ продолжай.")
        return 1

    print("[preflight] РЕЗУЛЬТАТ: ГОНКА ПРИ РЕКЛЕЙМЕ (3). Кто-то создал .lock в момент между "
          "удалением просроченной записи и нашей попыткой пересоздать её. Не трогать, "
          "сообщить и подождать.")
    return 3


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("role", help="ID роли, напр. A")
    parser.add_argument("--acquire", metavar="NICK", default=None,
                         help="попытаться атомарно занять лок под этим ником вместо read-only проверки")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT_MIN,
                         help=f"таймаут простоя в минутах (по умолчанию {DEFAULT_TIMEOUT_MIN})")
    args = parser.parse_args()

    role = args.role.upper()
    session_file = SESSIONS_DIR / f"{role}.md"
    if not session_file.exists():
        print(f"[preflight] ОШИБКА: {session_file} не найден (неизвестная роль?).")
        return 4

    if args.acquire:
        if not check_nicks_registered(args.acquire):
            print(f"[preflight] ПРЕДУПРЕЖДЕНИЕ: ник {args.acquire!r} не найден в NICKS.md "
                  "(зарегистрируй его ДО занятия лока по конвенции CHARTER §12, но это не блокирует .lock).")
        return acquire(role, args.acquire, args.timeout)

    return inspect(role, args.timeout)


if __name__ == "__main__":
    sys.exit(main())
