#!/usr/bin/env python3
"""Обновляет механические секции «Пульта оркестратора» (Claude Artifact) на месте.

Пульт — консоль в стиле IBM Plex ("что показывает git", "выработка по ролям",
"темп по неделям", "бюджет контекста", "журналы") + рукописные разделы
(срочное, что требуется от владельца, режимы работы) — вторые правит только
живая сессия ORCH, этот скрипт их не трогает.

Механические секции размечены HTML-комментариями `<!-- AUTO:<имя>:start -->`
/ `<!-- AUTO:<имя>:end -->` прямо в опубликованном артефакте. Скрипт читает
текущий HTML артефакта (передаётся файлом — Routine получает его через
Artifact({action:"read"}) и сохраняет на диск перед вызовом), вычисляет
свежие фрагменты из git и файлов репозитория, подменяет только то, что между
маркерами, и пишет результат — готовый к `Artifact` republish на тот же url.

Источники данных — те же, что у kpi_git.py (git log --all) и
check-context-budget.sh (coordination/roles/*.md против budget.json), функции
git-парсинга импортированы оттуда напрямую, а не переписаны заново.

Запуск (из корня репозитория):
    python coordination/tools/pult_update.py текущий_пульт.html --out новый_пульт.html
"""

import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpi_git import collect, aggregate, by_week  # переиспользуем парсинг git, не дублируем

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUDGET_JSON = os.path.join(ROOT, ".claude", "hooks", "budget.json")
ROLES_DIR = os.path.join(ROOT, "coordination", "roles")

# Исторические точки сравнения 15.08 — зафиксированы один раз, не пересчитываются
# (тот же принцип, что EXCLUDED_COMMITS в kpi_git.py: константа с объяснением,
# а не скрытое магическое число).
JOURNAL_BASELINES = {
    "ACTIVITY.md": 139_000,
    "HANDOFFS.md": 46_600,
    "QUESTIONS.md": 15_800,
}


def fmt(n):
    return f"{n:,}".replace(",", " ")


def ddmm(iso_date):
    """'2026-07-28' -> '28.07'."""
    y, m, d = iso_date.split("-")
    return f"{d}.{m}"


def fmt_k(tokens):
    """231000 → '231к'; 88100 → '88,1к' (запятая — русский десятичный разделитель)."""
    v = tokens / 1000
    if v >= 100:
        return f"{round(v)}к"
    return f"{v:.1f}".replace(".", ",") + "к"


def render_git_kpi(data):
    roles = data["roles"]
    named = sorted(
        ((r, v) for r, v in roles.items() if r not in ("(без тега)", "(merge)")),
        key=lambda kv: -kv[1]["added"],
    )
    untagged = roles.get("(без тега)")
    total_added = sum(v["added"] for v in roles.values())
    total_data = sum(v["data_added"] for v in roles.values())
    weeks = data.get("commits_by_week", {})

    tiles = f"""<div class="tiles">
  <div class="tile"><p class="k">коммитов учтено</p><p class="v">{data['commits_counted']}</p>
    <p class="n">два корня исключены — импорты дерева</p></div>
  <div class="tile"><p class="k">авторских строк</p><p class="v">{fmt(total_added)}</p>
    <p class="n">код, документы, записи</p></div>
  <div class="tile"><p class="k">внесено данных</p><p class="v">{fmt(total_data)}</p>
    <p class="n">спектры, вывод прогонов, архив</p></div>
  <div class="tile"><p class="k">активных недель</p><p class="v">{len(weeks)}</p>
    <p class="n">{ddmm(min(weeks)) if weeks else '?'} — {ddmm(max(weeks)) if weeks else '?'}</p></div>
</div>"""

    maxadd = named[0][1]["added"] if named else 1
    rows = []
    for role, v in named:
        pct = round(v["added"] / maxadd * 100, 1)
        rows.append(
            f"    <tr><td class='role'>{role}</td>"
            f"<td class='barcell'><div class='bar' style='width:{pct}%'></div></td>"
            f"<td class='num'>{fmt(v['added'])}</td><td class='num'>{v['commits']}</td>"
            f"<td class='num'>{v['active_days']}</td>"
            f"<td class='period'>{ddmm(v['first'])} – {ddmm(v['last'])}</td></tr>"
        )
    if untagged:
        rows.append(
            "    <tr><td class='role dim'>без тега</td>"
            "<td class='barcell'><span class='period'>не роль — вне сравнения</span></td>"
            f"<td class='num'>{fmt(untagged['added'])}</td><td class='num'>{untagged['commits']}</td>"
            f"<td class='num'>{untagged['active_days']}</td>"
            f"<td class='period'>{ddmm(untagged['first'])} – {ddmm(untagged['last'])}</td></tr>"
        )
    table = f"""<div class="scroll">
<table>
  <thead><tr>
    <th>роль</th><th>авторских строк</th><th class="num">строк</th><th class="num">комм.</th>
    <th class="num">дней</th><th>период</th>
  </tr></thead>
  <tbody>
{chr(10).join(rows)}
  </tbody>
</table>
</div>"""

    week_items = sorted(weeks.items())
    maxw = max(weeks.values()) if weeks else 1
    wk_divs = "\n    ".join(
        f"<div class='wk' title='неделя {ddmm(d)} — {n} коммит(ов)'><span class='lab'>{n}</span>"
        f"<div class='col' style='height:{round(n/maxw*100,1)}%'></div></div>"
        for d, n in week_items
    )
    wk_labels = "".join(f"<span>{ddmm(d)}</span>" for d, _ in week_items)
    chart = f"""<div class="chartbox">
  <div class="weeks">
    {wk_divs}
  </div>
  <div class="weeklabels">
    {wk_labels}
  </div>
</div>"""

    return tiles, table, chart


def render_budget():
    with open(BUDGET_JSON, "r", encoding="utf-8") as fh:
        budget = json.load(fh)
    limit_bytes = budget["files"][0]["limit_bytes"]

    sizes = {}
    for fn in sorted(os.listdir(ROLES_DIR)):
        if fn.endswith(".md"):
            sizes[fn] = os.path.getsize(os.path.join(ROLES_DIR, fn))
    if not sizes:
        return "<!-- нет roles/*.md -->"

    max_bytes = max(sizes.values())
    rows = []
    for fn, b in sorted(sizes.items(), key=lambda kv: -kv[1]):
        tok = round(b / 2)
        limit_tok = round(limit_bytes / 2)
        width = round(b / max_bytes * 100, 1)
        limit_pos = round(limit_bytes / max_bytes * 100, 1)
        ratio = round(b / limit_bytes, 1)
        over = b > limit_bytes
        cls = "ok" if not over else ("crit" if ratio >= 2.0 else "warn")
        fill_cls = "fill over" if over else "fill"
        rows.append(
            f"""  <div class="brow">
    <span class="fname">roles/{fn}</span>
    <div class="track"><div class="{fill_cls}" style="width:{width}%"></div>"""
            f"""<div class="limit" style="left:{limit_pos}%"></div></div>
    <span class="pill {cls}">×{ratio} · {fmt(tok)} / {fmt(limit_tok)}</span>
  </div>"""
        )
    return f'<div class="budget">\n{chr(10).join(rows)}\n</div>'


def render_journals():
    tiles = []
    for name, baseline in JOURNAL_BASELINES.items():
        path = os.path.join(ROOT, "coordination", name)
        b = os.path.getsize(path)
        tok = round(b / 2)
        ratio = round(tok / baseline, 1)
        tiles.append(
            f"  <div class='tile'><p class='k'>{name}</p><p class='v'>{fmt_k(tok)}</p>"
            f"<p class='n'>было {fmt_k(baseline)} · ×{str(ratio).replace('.', ',')}</p></div>"
        )
    tiles.append(
        "  <div class='tile'><p class='k'>база старта</p><p class='v'>33к</p>"
        "<p class='n'>неуправляема: промт и инструменты</p></div>"
    )
    return "<div class=\"tiles\">\n" + "\n".join(tiles) + "\n</div>"


def splice(html, name, fragment):
    pattern = re.compile(
        rf"(<!-- AUTO:{re.escape(name)}:start -->\n?).*?(\n?<!-- AUTO:{re.escape(name)}:end -->)",
        re.DOTALL,
    )
    new_html, n = pattern.subn(lambda m: m.group(1) + fragment + m.group(2), html)
    if n == 0:
        sys.exit(f"маркер AUTO:{name} не найден в HTML — секция не обновлена")
    return new_html


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("current_html", help="текущий HTML артефакта «Пульт оркестратора»")
    parser.add_argument("--out", required=True, help="куда записать обновлённый HTML")
    args = parser.parse_args()

    commits = collect()
    if not commits:
        sys.exit("git log --all не дал коммитов — запускать из корня репозитория")
    roles = aggregate(commits)
    weeks = by_week(commits)
    data = {
        "commits_counted": len(commits),
        "roles": roles,
        "commits_by_week": weeks,
    }

    with open(args.current_html, "r", encoding="utf-8") as fh:
        html = fh.read()

    tiles, table, chart = render_git_kpi(data)
    untagged_n = data["roles"].get("(без тега)", {}).get("commits", 0)
    git_kpi_fragment = (
        f"<section>\n<h2>Что показывает git</h2>\n"
        f"<p class=\"sub\">{data['commits_counted']} коммитов учтено (два корня исключены как "
        f"импорты дерева, не работа роли). Строки разделены на авторскую работу и внесённые "
        f"данные — без этого деления вклад роли A был бы завышен в разы.</p>\n\n{tiles}\n</section>\n\n"
        f"<section>\n<h2>Выработка по ролям</h2>\n"
        f"<p class=\"sub\">Столбик — авторские строки. Роли различаются подписью, а не цветом: "
        f"мера здесь одна.</p>\n\n{table}\n<p class=\"note\"><strong>«Без тега» — {untagged_n} "
        f"коммитов.</strong> Работа до конвенции префиксов и слияния веток. Приписать её ролям "
        f"нельзя, скрыть — значит потерять часть картины, поэтому строка стоит отдельно и "
        f"серым.</p>\n</section>\n\n"
        f"<section>\n<h2>Темп по неделям</h2>\n"
        f"<p class=\"sub\">Коммитов в неделю, понедельник — начало. Провал 09–12.08 — те коммиты "
        f"не сохранились нигде, они попали в импорт 13.08 как содержимое файлов, не как "
        f"история.</p>\n\n{chart}\n</section>"
    )
    html = splice(html, "git-kpi", git_kpi_fragment)
    html = splice(html, "budget", render_budget())
    html = splice(html, "journals", render_journals())

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"записано: {args.out}")


if __name__ == "__main__":
    main()
