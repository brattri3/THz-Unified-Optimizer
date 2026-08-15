#!/usr/bin/env python
"""Регрессионный стресс-тест атомарности preflight.py --acquire (CHARTER §12).

Запускает N НЕЗАВИСИМЫХ ОС-процессов почти одновременно, все бьются за один и тот же
coordination/sessions/<ROLE>.lock. Ожидание: ровно один получает код 0 (создал с нуля) или
1 (реклейм), остальные N-1 получают код 2 (занято) - НИ ОДНОГО lost update, НИ ОДНОГО кода 3.

Требует, чтобы .lock роли отсутствовал ДО запуска (роль FREE) - тест сам его не чистит после
себя, чтобы результат можно было проверить руками.

Usage: python coordination/tools/stress_test_lock.py [ROLE] [N]
"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PREFLIGHT = REPO_ROOT / "coordination" / "tools" / "preflight.py"
PY = sys.executable


def run_one(role: str, i: int):
    nick = f"stress-{i:03d}"
    proc = subprocess.run(
        [PY, str(PREFLIGHT), role, "--acquire", nick],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return nick, proc.returncode


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass

    role = sys.argv[1].upper() if len(sys.argv) > 1 else "A"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    lock_file = REPO_ROOT / "coordination" / "sessions" / f"{role}.lock"

    if lock_file.exists():
        print(f"ОШИБКА ПРЕДУСЛОВИЯ: {lock_file} уже существует, тест требует чистого FREE. Прерываю.")
        return 1

    print(f"Запускаю {n} параллельных процессов, все бьются за {lock_file} ...")
    results = []
    with ThreadPoolExecutor(max_workers=n) as pool:
        for nick, code in pool.map(lambda i: run_one(role, i), range(n)):
            results.append((nick, code))

    codes: dict[int, list[str]] = {}
    for nick, code in results:
        codes.setdefault(code, []).append(nick)

    print("\n=== Сводка по кодам возврата ===")
    for code in sorted(codes):
        holders = codes[code]
        preview = holders[:5] + (["..."] if len(holders) > 5 else [])
        print(f"  код {code}: {len(holders)} процессов -> {preview}")

    winners = codes.get(0, []) + codes.get(1, [])
    ok = True
    if len(winners) != 1:
        print(f"\nПРОВАЛ: победителей (код 0 или 1) должно быть РОВНО 1, а получилось {len(winners)}: {winners}")
        ok = False
    else:
        print(f"\nПобедитель: {winners[0]} (единственный, кто занял лок)")

    if 3 in codes:
        print(f"ПРОВАЛ: {len(codes[3])} процессов получили код 3 (неопределённое состояние) - не должно быть ни одного.")
        ok = False

    if lock_file.exists():
        holder = lock_file.read_text(encoding="utf-8").splitlines()[0].strip()
        print(f"\nСодержимое .lock сейчас: держатель = {holder!r}")
        if winners and holder != winners[0]:
            print(f"ПРОВАЛ: держатель в файле ({holder!r}) не совпадает с зафиксированным победителем ({winners[0]!r})!")
            ok = False
        elif winners and holder == winners[0]:
            print("OK: держатель в файле совпадает с единственным победителем.")

    print("\n=== ИТОГ:", "PASS" if ok else "FAIL", "===")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
