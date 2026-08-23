#!/usr/bin/env python3
"""KPI по git: выработка ролей, темп, распределение по времени.

Зачем именно git. Прежний сборщик (`kpi_stats.py`, теперь в archive/) читал локальные
транскрипты `~/.claude/projects/**/*.jsonl` и реестр `NICKS.md`. Обе опоры мертвы: реестр
закрыт 2026-08-07, транскрипты привязаны к конкретной машине и субагентов не видят вовсе.
git — единственный источник, который одинаково доступен с рабочего ПК, с ноутбука и из
облачного контейнера.

Чего этот сборщик НЕ делает и делать не может: он не знает стоимости в долларах и не знает
токенов. В git этих данных нет. Стоимость по ролям даёт только телеметрия
(CLAUDE_CODE_ENABLE_TELEMETRY + OTLP, метрика claude_code.cost.usage с атрибутом agent.name);
владелец решил 2026-08-23 её пока не поднимать. Здесь — выработка и темп, и ничего сверх.

Две ловушки, обе учтены (разбор — reports/GIT_HISTORY_2026-08-23.md):

1. У репозитория ДВА корня и нет общего предка. Считать по одному `main` бессмысленно —
   на нём не существует ролей B, D, L, P. Берём `--all`.
2. Оба корня — импорты рабочего дерева, а не работа роли: a1f69c2 (1243 файла, 618 872
   строки, помечен тегом [C]) и 91fc814 (250 файлов, 179 236 строк). Включить их — значит
   приписать одной роли всю историю проекта. Исключаются явно.

Запуск из корня репозитория:
    python coordination/tools/kpi_git.py                 # текстовая сводка
    python coordination/tools/kpi_git.py --json out.json # плюс машинный вывод
"""

import argparse
import collections
import datetime
import json
import re
import subprocess
import sys

# Импорты рабочего дерева, а не работа роли. См. шапку, ловушка 2.
EXCLUDED_COMMITS = {
    "a1f69c2537de0172a2d402dddf30cf6d308a7c90": "импорт при пересоздании репозитория 2026-08-13",
    "91fc814": "первичный импорт 2026-07-14",
}

ROLE_RE = re.compile(r"^\[([A-Z]+)\]")
MERGE_RE = re.compile(r"^Merge (branch|remote-tracking branch|pull request)")

# Пути, которые роль не пишет руками: измерения владельца, машинный вывод фитов, архив.
# Считать их выработкой — значит объявить рекордсменом того, кто закоммитил больше
# спектров. Проверено на этом репозитории: у роли A из 311 876 добавленных строк
# 176 311 — data_pool/*.txt (спектры по 651 строке), ещё 115 237 — research/results/**
# (JSON и PNG прогонов). Авторской работы остаётся ~20 тысяч, то есть метрика «строк»
# без этого деления завышала вклад A в пятнадцать раз.
NOT_AUTHORED = ("data_pool/", "archive/", "research/results/")


def is_authored(path):
    return not path.startswith(NOT_AUTHORED)

SEP = "\x1f"  # разделитель полей, в тексте коммита не встречается
REC = "\x1e"  # разделитель записей


def git(*args):
    """Вызов git; при ошибке — понятное сообщение, а не трейсбек."""
    try:
        out = subprocess.run(
            ("git",) + args, capture_output=True, text=True, check=True, encoding="utf-8"
        )
    except FileNotFoundError:
        sys.exit("git не найден в PATH")
    except subprocess.CalledProcessError as exc:
        sys.exit(f"git {' '.join(args)} → код {exc.returncode}\n{exc.stderr.strip()}")
    return out.stdout


def excluded(sha):
    """Полный sha или любой из настроенных префиксов."""
    return any(sha.startswith(key) for key in EXCLUDED_COMMITS)


def collect():
    """Один проход `git log --all --numstat` → записи по коммитам."""
    raw = git(
        "log", "--all", "--numstat",
        f"--format={REC}%H{SEP}%ad{SEP}%s", "--date=short",
    )

    commits = []
    for chunk in raw.split(REC):
        chunk = chunk.strip("\n")
        if not chunk:
            continue
        head, _, body = chunk.partition("\n")
        sha, date, subject = head.split(SEP, 2)
        if excluded(sha):
            continue

        added = removed = files = 0
        data_added = 0
        for line in body.splitlines():
            parts = line.split("\t")
            if len(parts) != 3:
                continue
            files += 1
            # "-" в numstat означает бинарный файл — строк у него нет
            a = int(parts[0]) if parts[0] != "-" else 0
            r = int(parts[1]) if parts[1] != "-" else 0
            if is_authored(parts[2]):
                added += a
                removed += r
            else:
                data_added += a

        match = ROLE_RE.match(subject)
        if match:
            role = match.group(1)
        elif MERGE_RE.match(subject):
            role = "(merge)"
        else:
            role = "(без тега)"

        commits.append(
            {
                "sha": sha[:7],
                "date": date,
                "role": role,
                "subject": subject,
                "added": added,
                "removed": removed,
                "data_added": data_added,
                "files": files,
            }
        )
    return commits


def aggregate(commits):
    roles = collections.defaultdict(
        lambda: {
            "commits": 0, "added": 0, "removed": 0, "data_added": 0,
            "files": 0, "days": set(),
        }
    )
    for c in commits:
        agg = roles[c["role"]]
        agg["commits"] += 1
        agg["added"] += c["added"]
        agg["removed"] += c["removed"]
        agg["data_added"] += c["data_added"]
        agg["files"] += c["files"]
        agg["days"].add(c["date"])

    out = {}
    for role, agg in roles.items():
        days = sorted(agg["days"])
        out[role] = {
            "commits": agg["commits"],
            "added": agg["added"],
            "removed": agg["removed"],
            "data_added": agg["data_added"],
            "files_touched": agg["files"],
            "active_days": len(days),
            "first": days[0],
            "last": days[-1],
            "lines_per_commit": round(agg["added"] / agg["commits"]) if agg["commits"] else 0,
        }
    return out


def by_week(commits):
    weeks = collections.Counter()
    for c in commits:
        day = datetime.date.fromisoformat(c["date"])
        monday = day - datetime.timedelta(days=day.weekday())
        weeks[monday.isoformat()] += 1
    return dict(sorted(weeks.items()))


def report(commits, roles, weeks):
    print("KPI по git — выработка ролей")
    print(f"собрано: {datetime.datetime.now():%Y-%m-%d %H:%M}")
    print()
    print(f"Коммитов учтено: {len(commits)}")
    excl = ", ".join(f"{sha[:7]} ({why})" for sha, why in EXCLUDED_COMMITS.items())
    print(f"Исключено как импорт: {excl}")
    print()

    head = (
        f"{'роль':<12}{'комм':>6}{'+работа':>10}{'-работа':>9}"
        f"{'+данные':>10}{'дней':>6}{'стр/комм':>10}  период"
    )
    print(head)
    print("-" * len(head))
    order = sorted(
        roles.items(),
        key=lambda kv: (kv[0].startswith("("), -kv[1]["commits"]),
    )
    for role, r in order:
        print(
            f"{role:<12}{r['commits']:>6}{r['added']:>10}{r['removed']:>9}"
            f"{r['data_added']:>10}{r['active_days']:>6}{r['lines_per_commit']:>10}"
            f"  {r['first']}..{r['last']}"
        )
    print()
    print("«работа» — код, документы, записи; «данные» — " + ", ".join(NOT_AUTHORED) + ":")
    print("измерения владельца, машинный вывод прогонов и архив роль вносит, а не пишет.")

    print()
    print("Коммитов по неделям (понедельник — начало недели):")
    peak = max(weeks.values()) if weeks else 1
    for monday, n in weeks.items():
        print(f"  {monday}  {'█' * max(1, round(n * 40 / peak)):<40} {n}")

    print()
    print("Чего здесь нет: стоимости и токенов — в git их не существует.")
    print("Стоимость по ролям даёт только телеметрия (claude_code.cost.usage, атрибут")
    print("agent.name). Разово посмотреть расход можно командами /usage и /insights.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", metavar="ФАЙЛ", help="записать машинный вывод")
    args = parser.parse_args()

    commits = collect()
    if not commits:
        sys.exit("ни одного коммита не учтено — запускать из корня репозитория")

    roles = aggregate(commits)
    weeks = by_week(commits)
    report(commits, roles, weeks)

    if args.json:
        payload = {
            "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "excluded_commits": EXCLUDED_COMMITS,
            "commits_counted": len(commits),
            "roles": roles,
            "commits_by_week": weeks,
            "caveat": "git не содержит стоимости и токенов; выработка и темп — всё, что здесь есть",
        }
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2)
        print(f"\nмашинный вывод: {args.json}")


if __name__ == "__main__":
    main()
