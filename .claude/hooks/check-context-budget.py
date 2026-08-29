#!/usr/bin/env python3
"""Хук бюджета контекста холодного старта (PARADIGM_REVIEW_2026-08-23 §7 п.1).

Зачем: калибровка "цена в CLAUDE.md" протухала молча (Д-2) — файлы росли, таблица нет, роль
решала "читать/не читать" по неверной цифре. Правило "ужать до 1200 токенов" из
COLDSTART_REDESIGN_2026-08-15 не проверялось нигде. Правила, встроенные в действие, держатся
сами (CHARTER §4) — этот хук переводит "не забыть" в "проверяется", как уже сделано для
трейлеров. Только предупреждает, никогда не блокирует.

Единица — байты UTF-8 (`wc -c`), не Unicode code points и не токены: это то, чем на самом деле
измерена вся история проекта (roles/A.md — 9739 байт ≈ 4870 ток. × 2.0 симв/ток из калибровки —
сходится с байтами, не с числом символов). См. budget.json.

⚠ Меряем в LF, а не в том, что лежит на диске (правка 26.08). Файлы ужимались в облачном контуре
(LF), а на рабочем ПК владельца git раскладывает их с CRLF — лишний байт на строку, ~35 на файл
роли. Из-за этого шесть из семи roles/*.md, ужатые ровно под лимит, локально показывали
2418–2437 и все считались нарушителями. Лимит был верен, врало измерение; предупреждение из-за
этого читалось как шум, а хук ровно за тем и заведён, чтобы ему верили. CR отбрасываем перед
счётом.

27.08: перенесено с bash на Python при сверке с последней версией multi-agent-coordination-skill
(тот же скилл, извлечённый из этого же проекта 24.08, с тех пор доработан отдельно, ветка
fix/issue-6-dashboard-integrity — на момент переноса ещё не смержена в её собственный main).
Перенос закрывает реальный пробел bash-версии, не найденный в CRLF-инциденте: `budget.json`
там был декоративен — хук вытаскивал `limit_bytes` из ПЕРВОГО совпадения через grep и молча
гонял только захардкоженный `coordination/roles/*.md`, игнорируя остальные записи `files`, если
бы они появились. Сейчас в budget.json ровно одно правило, поэтому баг был незаметен — но
`files` уже array, рассчитан на несколько правил, и молчаливо игнорировать вторую запись —
ровно то, от чего этот хук должен защищать. Поддержка `git worktree` (`.git` — файл, не каталог)
у bash-версии уже была корректной через `git rev-parse --show-toplevel` — переносить было нечего.

Событие: SessionStart, только source=startup|clear — resume/compact файлы старта заново не
читают, предупреждать там — чистый шум.
"""

import argparse
import fnmatch
import json
import os
import subprocess
import sys

DEFAULT_LIMIT = 2400
DEFAULT_GLOB = "coordination/roles/*.md"


def get_project_root(start_dir):
    """Подняться до корня проекта.

    Сначала git, потом каталог с .claude, потом каталог с .git.

    Проверка `.git` — через exists(), не isdir(): внутри `git worktree` (CHARTER §5 рекомендует
    его для активно правящих общий код) `.git` — ФАЙЛ, а не каталог.
    """
    start = os.path.abspath(start_dir)
    try:
        result = subprocess.run(
            ["git", "-C", start, "rev-parse", "--show-toplevel"],
            capture_output=True, text=True, check=True,
        )
        root = result.stdout.strip()
        if root and os.path.isdir(root):
            return root
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        pass

    current = start
    while True:
        if os.path.isdir(os.path.join(current, ".claude")):
            return current
        if os.path.exists(os.path.join(current, ".git")):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            return start
        current = parent


def load_rules(root):
    """Прочитать budget.json в [(glob, limit_bytes), ...].

    Каждая запись `files` учитывается по-настоящему, и glob используется КАК glob — не как
    источник одного числа для захардкоженного пути (см. историю правки в шапке файла).
    """
    config_path = os.path.join(root, ".claude", "hooks", "budget.json")
    if not os.path.exists(config_path):
        return [(DEFAULT_GLOB, DEFAULT_LIMIT)]
    try:
        with open(config_path, "r", encoding="utf-8") as handle:
            cfg = json.load(handle)
    except (OSError, ValueError):
        return [(DEFAULT_GLOB, DEFAULT_LIMIT)]

    rules = []
    for rule in cfg.get("files", []):
        pattern = rule.get("glob")
        if not pattern:
            continue
        rules.append((pattern, rule.get("limit_bytes", DEFAULT_LIMIT)))
    return rules or [(DEFAULT_GLOB, DEFAULT_LIMIT)]


def iter_matches(root, pattern):
    """Файлы под `root`, совпавшие с POSIX-glob (включая `**`)."""
    normalised = pattern.replace("\\", "/")
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in (".git", "node_modules", "__pycache__")]
        for name in filenames:
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, root).replace(os.sep, "/")
            if fnmatch.fnmatch(rel, normalised):
                yield rel, full


def measure(path):
    """Размер в байтах с нормализацией CRLF→LF, или None, если файл не читается."""
    try:
        with open(path, "rb") as handle:
            return len(handle.read().replace(b"\r\n", b"\n"))
    except OSError:
        return None


def should_run(argv_forced):
    """True, когда именно этот запуск должен проверять.

    SessionStart срабатывает на несколько источников; только startup и clear заново читают файл
    холодного старта — предупреждать на resume/compact — чистый шум. Всё неразбираемое
    запускает проверку: для предупреждения безопаснее упасть в сторону "проверить".
    """
    if argv_forced:
        return True
    try:
        if sys.stdin is None or sys.stdin.isatty():
            return True
        payload = json.loads(sys.stdin.read())
    except Exception:
        return True
    return payload.get("source") in ("startup", "clear")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strict", action="store_true",
                     help="код возврата 2 при превышении (по умолчанию: только предупреждение, 0)")
    ap.add_argument("--json", action="store_true",
                     help="вывести находки как JSON в stdout (для CI)")
    ap.add_argument("--force", action="store_true",
                     help="проверить независимо от source в stdin SessionStart")
    args = ap.parse_args()

    if not should_run(args.force):
        return 0

    root = get_project_root(os.getcwd())
    oversized = []
    unreadable = []

    for pattern, limit in load_rules(root):
        for rel, full in iter_matches(root, pattern):
            size = measure(full)
            if size is None:
                unreadable.append(rel)
            elif size > limit:
                oversized.append({"file": rel, "bytes": size, "limit": limit})

    if args.json:
        print(json.dumps({"oversized": oversized, "unreadable": unreadable}, ensure_ascii=False))
        return 2 if (args.strict and oversized) else 0

    if not oversized and not unreadable:
        return 0

    parts = []
    if oversized:
        listed = ", ".join(
            f"{item['file']} ({item['bytes']}/{item['limit']} Б)" for item in oversized
        )
        parts.append(f"Бюджет холодного старта превышен: {listed}.")
    if unreadable:
        parts.append("Не удалось прочитать: " + ", ".join(unreadable) + ".")

    parts.append("Лимит coordination/roles/*.md — из budget.json.")

    payload = json.dumps({"systemMessage": " ".join(parts)}, ensure_ascii=False)
    sys.stdout.buffer.write(payload.encode("utf-8") + b"\n")
    sys.stdout.buffer.flush()

    return 2 if (args.strict and oversized) else 0


if __name__ == "__main__":
    raise SystemExit(main())
