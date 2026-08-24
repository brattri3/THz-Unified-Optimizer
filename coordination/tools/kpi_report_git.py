#!/usr/bin/env python3
"""Рендерер HTML-отчёта по выработке ролей из вывода kpi_git.py.

Берёт JSON (`kpi_git.py --json ...`) и строит самодостаточный HTML: карточки по
ролям, три графика (коммиты по ролям, авторская работа по ролям, коммиты по
неделям), полная таблица, методология. Без внешних зависимостей — только stdlib,
без CDN, тема light/dark берётся из `prefers-color-scheme` + `[data-theme]`.

Категориальная палитра — 7 первых слотов валидированной темы дизайн-скилла
dataviz (см. references/palette.md), проверена валидатором на поверхностях
этого отчёта (light #ffffff, dark #171a1e) — оба режима проходят все проверки.
Роли вне этого списка (новая роль, ранее не встречавшаяся) получают приглушённый
резервный цвет — это осознанный выбор дизайна, не баг: девятая и далее серия не
получает выдуманный оттенок (см. dataviz/anti-patterns).

Запуск:
    python coordination/tools/kpi_git.py --json /tmp/kpi_git.json
    python coordination/tools/kpi_report_git.py /tmp/kpi_git.json --out coordination/reports/kpi_report_git.html
"""

import argparse
import html
import json
import sys

# Роль → (light, dark) — первые 7 слотов валидированной категориальной темы
# dataviz (blue/orange/aqua/yellow/magenta/green/violet), проверено
# validate_palette.js на поверхностях этого отчёта. Порядок фиксирован по
# идентичности роли, не по рангу выработки (CVD-безопасность зависит от
# порядка слотов, см. references/palette.md).
ROLE_COLOR = {
    "ORCH": ("#2a78d6", "#3987e5"),
    "A":    ("#eb6834", "#d95926"),
    "B":    ("#1baf7a", "#199e70"),
    "C":    ("#eda100", "#c98500"),
    "D":    ("#e87ba4", "#d55181"),
    "L":    ("#008300", "#008300"),
    "P":    ("#4a3aa7", "#9085e9"),
}
# Резерв для роли, которой ещё нет в ROLE_COLOR (новая роль в проекте) —
# приглушённый, а не выдуманный яркий оттенок.
FALLBACK_COLOR = ("#7a756a", "#9b968a")
# "(без тега)"/"(merge)" — не роли, приглушённый нейтральный цвет намеренно
# (dataviz: девятая+ серия не получает свой оттенок, сворачивается в "прочее").
OTHER_COLOR = ("#a39d8e", "#6d6656")


def esc(s):
    return html.escape(str(s), quote=True)


def role_color(role, dark):
    if role in ("(без тега)", "(merge)"):
        return OTHER_COLOR[1 if dark else 0]
    return ROLE_COLOR.get(role, FALLBACK_COLOR)[1 if dark else 0]


def fmt(n):
    return f"{n:,}".replace(",", " ")


def bar_chart(items, value_key, title, unit=""):
    """items: [(role, {..данные..})], value_key: имя числового поля."""
    if not items:
        return "<p class='subtitle'>Нет данных.</p>"
    maxv = max(it[1][value_key] for it in items) or 1
    row_h, gap, left_pad, right_pad = 26, 12, 8, 56
    width = 640
    bar_max_w = width - left_pad - right_pad - 90  # 90 — под подпись роли слева
    height = len(items) * (row_h + gap) + 8
    rows = []
    for i, (role, data) in enumerate(items):
        v = data[value_key]
        y = i * (row_h + gap) + 4
        w = max(2, round(v / maxv * bar_max_w))
        rows.append(
            f"<text x='0' y='{y + row_h * 0.68:.1f}' class='chart-label' "
            f"font-weight='700'>{esc(role)}</text>"
            f"<rect x='90' y='{y}' width='{w}' height='{row_h - 4}' rx='4' "
            f"fill='var(--role-{esc(role_key(role))})'>"
            f"<title>{esc(role)}: {fmt(v)}{esc(unit)}</title></rect>"
            f"<text x='{90 + w + 8}' y='{y + row_h * 0.68:.1f}' "
            f"class='chart-label tnum'>{fmt(v)}{esc(unit)}</text>"
        )
    return (
        f"<svg viewBox='0 0 {width} {height}' width='100%' height='{height}' "
        f"class='chart' role='img' aria-label='{esc(title)}'>{''.join(rows)}</svg>"
    )


def role_key(role):
    """Ключ CSS-переменной — только буквы, спецкатегории → 'other'."""
    return role if role.isalpha() and role.isupper() else "other"


def week_chart(weeks):
    if not weeks:
        return "<p class='subtitle'>Нет данных.</p>"
    items = sorted(weeks.items())
    maxv = max(v for _, v in items) or 1
    width, height = 640, 220
    left_pad, bottom_pad, top_pad = 34, 28, 10
    plot_h = height - top_pad - bottom_pad
    bar_w = (width - left_pad - 8) / len(items)
    rows = [
        f"<line x1='{left_pad}' y1='{top_pad}' x2='{left_pad}' "
        f"y2='{top_pad + plot_h}' class='grid-line'/>",
        f"<line x1='{left_pad}' y1='{top_pad + plot_h}' x2='{width}' "
        f"y2='{top_pad + plot_h}' class='grid-line'/>",
    ]
    for i, (monday, n) in enumerate(items):
        bh = round(n / maxv * (plot_h - 4))
        x = left_pad + i * bar_w + 3
        y = top_pad + plot_h - bh
        rows.append(
            f"<rect x='{x:.1f}' y='{y}' width='{max(2, bar_w - 6):.1f}' "
            f"height='{bh}' rx='3' fill='var(--accent)'>"
            f"<title>неделя {esc(monday)}: {n} коммит(ов)</title></rect>"
        )
        if len(items) <= 12 or i % max(1, len(items) // 10) == 0:
            rows.append(
                f"<text x='{x + bar_w / 2:.1f}' y='{height - 6}' "
                f"class='chart-legend' text-anchor='middle'>{esc(monday[5:])}</text>"
            )
    rows.append(
        f"<text x='4' y='{top_pad + 8}' class='chart-legend'>{maxv}</text>"
    )
    return (
        f"<svg viewBox='0 0 {width} {height}' width='100%' height='{height}' "
        f"class='chart' role='img' aria-label='Коммиты по неделям'>{''.join(rows)}</svg>"
    )


def role_vars(dark):
    out = []
    for role in ROLE_COLOR:
        out.append(f"--role-{role}:{role_color(role, dark)};")
    out.append(f"--role-other:{OTHER_COLOR[1 if dark else 0]};")
    return "".join(out)


CSS = f"""
:root{{
  --bg:#f7f4ee; --panel:#ffffff; --border:#e4ded2; --fg:#1c1b18; --muted:#6d6656;
  --accent:#b8650f; --code-bg:#efe9db; --grid-line:#e4ded2;
  --callout-bg:#fdf1dc; --callout-fg:#8a5410; --callout-border:#eacea0;
  {role_vars(False)}
}}
@media (prefers-color-scheme: dark){{
  :root:not([data-theme="light"]){{
    --bg:#0f1113; --panel:#171a1e; --border:#2a2e33; --fg:#e9ece7; --muted:#8b9199;
    --accent:#ffb545; --code-bg:#1d2126; --grid-line:#262a2f;
    --callout-bg:#241c0c; --callout-fg:#f0c46b; --callout-border:#4a3a17;
    {role_vars(True)}
  }}
}}
:root[data-theme="dark"]{{
  --bg:#0f1113; --panel:#171a1e; --border:#2a2e33; --fg:#e9ece7; --muted:#8b9199;
  --accent:#ffb545; --code-bg:#1d2126; --grid-line:#262a2f;
  --callout-bg:#241c0c; --callout-fg:#f0c46b; --callout-border:#4a3a17;
  {role_vars(True)}
}}
:root[data-theme="light"]{{
  --bg:#f7f4ee; --panel:#ffffff; --border:#e4ded2; --fg:#1c1b18; --muted:#6d6656;
  --accent:#b8650f; --code-bg:#efe9db; --grid-line:#e4ded2;
  --callout-bg:#fdf1dc; --callout-fg:#8a5410; --callout-border:#eacea0;
  {role_vars(False)}
}}
*{{box-sizing:border-box;}}
body{{background:var(--bg); color:var(--fg); margin:0;}}
.tnum{{font-variant-numeric:tabular-nums;}}
.wrap{{
  max-width:960px; margin:0 auto; padding:28px 18px 72px;
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif; line-height:1.55;
}}
h1{{
  font-family:ui-monospace,SFMono-Regular,"Cascadia Code",Consolas,monospace;
  font-size:1.4rem; font-weight:700; margin:0 0 4px; text-wrap:balance;
}}
h2{{
  font-family:ui-monospace,SFMono-Regular,"Cascadia Code",Consolas,monospace;
  font-size:.76rem; font-weight:700; text-transform:uppercase; letter-spacing:.08em;
  color:var(--muted); margin:0 0 12px;
}}
section{{margin-top:32px;}}
.subtitle{{color:var(--muted); font-size:.86rem; margin:0 0 4px;}}
.callout{{
  background:var(--callout-bg); color:var(--callout-fg); border:1px solid var(--callout-border);
  border-radius:6px; padding:11px 15px; font-size:.85rem; margin:14px 0 0;
}}
.callout code{{background:transparent; padding:0; color:inherit; border-bottom:1px dotted currentColor;}}
.cards{{display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:12px;}}
.card{{
  background:var(--panel); border:1px solid var(--border); border-top:3px solid var(--role-color);
  border-radius:8px; padding:13px 14px 11px;
}}
.card .role{{font-weight:700; font-family:ui-monospace,monospace; font-size:.9rem;}}
.card .metric{{font-size:1.28rem; font-weight:700; font-family:ui-monospace,monospace; margin-top:6px;}}
.card .label{{font-size:.66rem; color:var(--muted); text-transform:uppercase; letter-spacing:.04em;}}
.card .period{{font-size:.72rem; color:var(--muted); margin-top:8px;}}
.panel{{background:var(--panel); border:1px solid var(--border); border-radius:8px; padding:16px 18px 14px;}}
.chart-label{{fill:var(--fg); font-size:11.5px; font-family:ui-monospace,monospace;}}
.chart-legend{{fill:var(--muted); font-size:10px; font-family:ui-monospace,monospace;}}
.grid-line{{stroke:var(--grid-line); stroke-width:1;}}
.table-scroll{{overflow-x:auto; border:1px solid var(--border); border-radius:8px;}}
.data-table{{width:100%; border-collapse:collapse; font-size:.82rem; min-width:560px;}}
.data-table th{{
  text-align:left; color:var(--muted); font-weight:600; padding:8px 10px;
  border-bottom:1px solid var(--border); font-size:.66rem; text-transform:uppercase;
  letter-spacing:.04em; white-space:nowrap;
}}
.data-table td{{padding:7px 10px; border-bottom:1px solid var(--border); font-variant-numeric:tabular-nums;}}
.data-table tbody tr:hover{{background:color-mix(in srgb, var(--accent) 6%, transparent);}}
.role-cell{{border-left:4px solid var(--role-color); font-weight:700; font-family:ui-monospace,monospace;}}
code{{background:var(--code-bg); padding:1px 6px; border-radius:4px; font-size:.85em; font-family:ui-monospace,monospace;}}
ul.method{{padding-left:1.1em; display:flex; flex-direction:column; gap:7px; font-size:.86rem;}}
footer{{color:var(--muted); font-size:.75rem; margin-top:40px; border-top:1px solid var(--border); padding-top:12px;}}
a{{color:var(--accent);}}
:focus-visible{{outline:2px solid var(--accent); outline-offset:2px;}}
"""


def build_html(data):
    roles = data["roles"]
    # Именованные роли — карточки и графики; "(без тега)"/"(merge)" — только таблица.
    named = sorted(
        ((r, v) for r, v in roles.items() if r not in ("(без тега)", "(merge)")),
        key=lambda kv: -kv[1]["commits"],
    )
    all_rows = sorted(
        roles.items(),
        key=lambda kv: (kv[0] in ("(без тега)", "(merge)"), -kv[1]["commits"]),
    )

    cards = "".join(
        f"<div class='card' style='--role-color:var(--role-{esc(role_key(r))})'>"
        f"<div class='role'>{esc(r)}</div>"
        f"<div class='metric tnum'>{v['commits']}</div><div class='label'>коммитов</div>"
        f"<div class='metric tnum' style='font-size:1rem;margin-top:8px'>+{fmt(v['added'])}"
        f"/−{fmt(v['removed'])}</div><div class='label'>строк работы</div>"
        f"<div class='period'>{esc(v['first'])} … {esc(v['last'])} · {v['active_days']} дн.</div>"
        f"</div>"
        for r, v in named
    )

    commits_chart = bar_chart(named, "commits", "Коммиты по ролям")
    work_chart = bar_chart(
        sorted(named, key=lambda kv: -kv[1]["added"]), "added", "Авторская работа по ролям", " стр."
    )
    weeks_chart = week_chart(data.get("commits_by_week", {}))

    table_rows = "".join(
        f"<tr><td class='role-cell' style='--role-color:var(--role-{esc(role_key(r))})'>{esc(r)}</td>"
        f"<td class='tnum'>{v['commits']}</td>"
        f"<td class='tnum'>+{fmt(v['added'])}</td><td class='tnum'>−{fmt(v['removed'])}</td>"
        f"<td class='tnum'>{fmt(v['data_added'])}</td>"
        f"<td class='tnum'>{v['active_days']}</td>"
        f"<td class='tnum'>{v['lines_per_commit']}</td>"
        f"<td>{esc(v['first'])} … {esc(v['last'])}</td></tr>"
        for r, v in all_rows
    )

    excluded = "; ".join(
        f"<code>{esc(sha[:7])}</code> ({esc(why)})"
        for sha, why in data.get("excluded_commits", {}).items()
    )

    return f"""<title>THz KPI по git</title>
<style>{CSS}</style>
<div class="wrap">
  <h1>Выработка ролей &middot; THz-Unified-Optimizer</h1>
  <div class="subtitle">сгенерировано {esc(data.get('generated_at', '?'))} &middot;
    источник &mdash; <code>git log --all --numstat</code>, единственный, одинаково доступный
    во всех трёх контурах (рабочий ПК, ноутбук, облако)</div>
  <div class="callout">
    ⚠ Здесь нет стоимости и токенов &mdash; в git их не существует. Стоимость по ролям даёт только
    телеметрия (<code>claude_code.cost.usage</code>, атрибут <code>agent.name</code>) &mdash; не
    поднята (см. <code>PARADIGM_REVIEW_2026-08-23.md</code> Д-5). Разовый срез расхода &mdash;
    команды <code>/usage</code> и <code>/insights</code>.
  </div>

  <section>
    <h2>Сводка по ролям</h2>
    <div class="cards">{cards}</div>
  </section>

  <section>
    <h2>Коммиты по ролям</h2>
    <div class="panel">{commits_chart}</div>
  </section>

  <section>
    <h2>Авторская работа (строк) по ролям</h2>
    <p class="subtitle">Без данных измерений/машинного вывода &mdash; см. методологию.</p>
    <div class="panel">{work_chart}</div>
  </section>

  <section>
    <h2>Коммиты по неделям</h2>
    <div class="panel">{weeks_chart}</div>
  </section>

  <section>
    <h2>Полная таблица</h2>
    <div class="table-scroll"><table class="data-table">
      <thead><tr><th>Роль</th><th>Комм.</th><th>+работа</th><th>−работа</th><th>+данные</th>
        <th>Дней</th><th>Стр/комм</th><th>Период</th></tr></thead>
      <tbody>{table_rows}</tbody>
    </table></div>
  </section>

  <section>
    <h2>Методология и оговорки</h2>
    <ul class="method">
      <li><b>Считает по <code>--all</code></b>, не по <code>main</code> &mdash; у репозитория два
        корня без общего предка, на <code>main</code> не существует части ролей.</li>
      <li><b>Исключены корневые импорты</b> рабочего дерева, не работа роли: {excluded or '&mdash;'}.</li>
      <li><b>«Работа» vs «данные»</b> &mdash; строки в <code>data_pool/</code>, <code>archive/</code>,
        <code>research/results/</code> считаются отдельно от авторской работы: без этого деления
        роль, коммитящая много измерений/прогонов, выглядит кратно продуктивнее.</li>
      <li><b>«(без тега)»/«(merge)»</b> &mdash; коммиты без префикса <code>[РОЛЬ]</code> и мержи;
        не роли, в графики и карточки не входят, только в полную таблицу.</li>
      <li><b>Чего здесь нет:</b> стоимости и токенов (не в git), статуса «кто сейчас работает»
        (на это отвечает <code>claude agents</code>, отдельный реестр процессов этот проект уже
        один раз отменил &mdash; В-27).</li>
    </ul>
  </section>

  <footer>
    Сгенерировано <code>coordination/tools/kpi_report_git.py</code> из вывода
    <code>coordination/tools/kpi_git.py --json</code>. Обновляется по расписанию (Routine).
    Ручной перезапуск &mdash; из корня репозитория:
    <code>python coordination/tools/kpi_git.py --json /tmp/k.json &amp;&amp;
    python coordination/tools/kpi_report_git.py /tmp/k.json --out coordination/reports/kpi_report_git.html</code>
  </footer>
</div>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("json_in", help="вывод kpi_git.py --json")
    parser.add_argument("--out", required=True, help="куда записать HTML")
    args = parser.parse_args()

    with open(args.json_in, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not data.get("roles"):
        sys.exit("во входном JSON нет данных по ролям")

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(build_html(data))
    print(f"записано: {args.out}")


if __name__ == "__main__":
    main()
