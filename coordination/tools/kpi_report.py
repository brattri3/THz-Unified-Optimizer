"""
Рендер HTML-отчёта из данных `kpi_stats.build_data()`. Никаких внешних библиотек и CDN -
графики - это SVG, сгенерированный прямо здесь строками, чтобы отчёт был самодостаточным
файлом (можно открыть локально в браузере или опубликовать как Artifact).

Вызывается из kpi_stats.py: `kpi_report.render(data) -> str` (HTML-фрагмент, без
<!doctype>/<html>/<head>/<body> - их добавляет либо браузер неявно, либо Artifact-обёртка).
"""
from __future__ import annotations

import html
from datetime import datetime


def _fmt_int(n: float) -> str:
    return f"{int(round(n)):,}".replace(",", " ")


def _fmt_money(x: float) -> str:
    return f"${x:,.2f}".replace(",", " ")


def _fmt_h(h: float) -> str:
    if h < 48:
        return f"{h:.1f} ч"
    return f"{h/24:.1f} дн"


def _esc(s: str) -> str:
    return html.escape(str(s))


# ---------------------------------------------------------------------------
# SVG-примитивы
# ---------------------------------------------------------------------------

def svg_stacked_tokens(roles: dict, role_order: list, role_color: dict, width=720) -> str:
    """Горизонтальные стековые бары: input/output/cache_write/cache_read по ролям."""
    seg_colors = {
        "input": "#6d7280", "output": "#d98a3d",
        "cache_write": "#8a6b9e", "cache_read": "#4f8fae",
    }
    seg_labels = {"input": "input", "output": "output", "cache_write": "cache-write", "cache_read": "cache-read"}
    rows = [(r, roles[r]) for r in role_order if roles[r]["tokens_grand_total"] > 0]
    if not rows:
        return "<p><em>Нет данных о токенах.</em></p>"
    max_total = max(r["tokens_grand_total"] for _, r in rows) or 1
    bar_h = 26
    gap = 14
    left_pad = 60
    right_pad = 170
    plot_w = width - left_pad - right_pad
    top_pad = 22
    height = len(rows) * (bar_h + gap) + gap + top_pad + 30
    svg = [f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" class="chart" role="img" '
           f'aria-label="Токены по ролям">']
    plot_bottom = top_pad + len(rows) * (bar_h + gap)
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        gx = left_pad + frac * plot_w
        svg.append(f'<line x1="{gx:.1f}" y1="{top_pad-6}" x2="{gx:.1f}" y2="{plot_bottom}" class="grid-line"/>')
        svg.append(f'<text x="{gx:.1f}" y="12" class="chart-legend" text-anchor="middle">{_fmt_int(max_total*frac)}</text>')
    y = top_pad
    for role, r in rows:
        x = left_pad
        for key in ("input", "output", "cache_write", "cache_read"):
            val = r["tokens_total"][key]
            w = (val / max_total) * plot_w
            if w > 0:
                svg.append(
                    f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{seg_colors[key]}">'
                    f'<title>{role} · {seg_labels[key]}: {_fmt_int(val)} ток.</title></rect>'
                )
                x += w
        svg.append(f'<text x="8" y="{y + bar_h/2 + 5}" class="chart-label" font-weight="700">{_esc(role)}</text>')
        svg.append(f'<text x="{left_pad + plot_w + 10}" y="{y + bar_h/2 + 5}" class="chart-label">'
                    f'{_fmt_int(r["tokens_grand_total"])}</text>')
        y += bar_h + gap
    # легенда
    lx = left_pad
    ly = height - 8
    for key in ("input", "output", "cache_write", "cache_read"):
        svg.append(f'<rect x="{lx}" y="{ly-10}" width="12" height="12" fill="{seg_colors[key]}"/>')
        svg.append(f'<text x="{lx+16}" y="{ly}" class="chart-legend">{seg_labels[key]}</text>')
        lx += 110
    svg.append("</svg>")
    return "".join(svg)


def svg_simple_bars(items: list[tuple[str, float]], role_color: dict, width=720, fmt=_fmt_int, unit="") -> str:
    if not items:
        return "<p><em>Нет данных.</em></p>"
    max_v = max(v for _, v in items) or 1
    bar_h = 26
    gap = 14
    left_pad = 60
    right_pad = 130
    plot_w = width - left_pad - right_pad
    top_pad = 22
    height = len(items) * (bar_h + gap) + gap + top_pad
    svg = [f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" class="chart" role="img">']
    plot_bottom = top_pad + len(items) * (bar_h + gap) - gap + bar_h
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        gx = left_pad + frac * plot_w
        svg.append(f'<line x1="{gx:.1f}" y1="{top_pad-6}" x2="{gx:.1f}" y2="{plot_bottom}" class="grid-line"/>')
        svg.append(f'<text x="{gx:.1f}" y="12" class="chart-legend" text-anchor="middle">{fmt(max_v*frac)}</text>')
    y = top_pad
    for role, val in items:
        w = (val / max_v) * plot_w if max_v else 0
        color = role_color.get(role, "#888")
        svg.append(f'<rect x="{left_pad}" y="{y}" width="{w:.1f}" height="{bar_h}" fill="{color}" rx="3">'
                    f'<title>{role}: {fmt(val)}{unit}</title></rect>')
        svg.append(f'<text x="8" y="{y + bar_h/2 + 5}" class="chart-label" font-weight="700">{_esc(role)}</text>')
        svg.append(f'<text x="{left_pad + w + 8}" y="{y + bar_h/2 + 5}" class="chart-label">{fmt(val)}{unit}</text>')
        y += bar_h + gap
    svg.append("</svg>")
    return "".join(svg)


def svg_gantt(chats: list[dict], role_color: dict, width=760) -> str:
    found = [c for c in chats if c.get("found") and c.get("first_ts") and c.get("last_ts")]
    if not found:
        return "<p><em>Нет данных о времени жизни чатов.</em></p>"
    parsed = []
    for c in found:
        f = datetime.fromisoformat(c["first_ts"])
        l = datetime.fromisoformat(c["last_ts"])
        parsed.append((c, f, l))
    parsed.sort(key=lambda x: x[1])
    t_min = min(p[1] for p in parsed)
    t_max = max(p[2] for p in parsed)
    span = (t_max - t_min).total_seconds() or 1

    row_h = 24
    gap = 6
    left_pad = 70
    right_pad = 20
    top_pad = 10
    plot_w = width - left_pad - right_pad
    height = len(parsed) * (row_h + gap) + top_pad + 30

    svg = [f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" class="chart" role="img" '
           f'aria-label="Хронология жизни чатов">']
    # сетка по дням
    import datetime as _dt
    day = t_min.replace(hour=0, minute=0, second=0, microsecond=0)
    while day <= t_max:
        x = left_pad + ((day - t_min).total_seconds() / span) * plot_w
        svg.append(f'<line x1="{x:.1f}" y1="{top_pad}" x2="{x:.1f}" y2="{height-22}" class="grid-line"/>')
        svg.append(f'<text x="{x:.1f}" y="{height-8}" class="chart-legend" text-anchor="middle">{day.strftime("%d.%m")}</text>')
        day += _dt.timedelta(days=1)

    y = top_pad
    for c, f, l in parsed:
        x = left_pad + ((f - t_min).total_seconds() / span) * plot_w
        w = max(((l - f).total_seconds() / span) * plot_w, 3)
        color = role_color.get(c["role"], "#888")
        dur = (l - f).total_seconds() / 3600
        svg.append(
            f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{row_h}" fill="{color}" rx="4">'
            f'<title>{c["role"]}/{c["nick"]}: {f.strftime("%Y-%m-%d %H:%M")} → '
            f'{l.strftime("%Y-%m-%d %H:%M")} ({dur:.1f} ч)</title></rect>'
        )
        svg.append(f'<text x="4" y="{y + row_h/2 + 4}" class="chart-label">{_esc(c["nick"])}</text>')
        y += row_h + gap
    svg.append("</svg>")
    return "".join(svg)


# ---------------------------------------------------------------------------
# Таблицы
# ---------------------------------------------------------------------------

def table_efficiency(data: dict) -> str:
    rows = []
    for r in data["role_order"]:
        role = data["roles"][r]
        git = data["git"].get(r, {"commits": 0, "insertions": 0, "deletions": 0})
        rlog = data["research_log"].get(r, 0)
        tokens = role["tokens_grand_total"]
        commits = git["commits"]
        tok_per_commit = tokens / commits if commits else None
        cost_per_rlog = role["cost_usd"] / rlog if rlog else None
        tool_per_turn = role["tool_calls_total"] / role["user_turns"] if role["user_turns"] else None
        rows.append((r, tokens, commits, git["insertions"], git["deletions"], rlog,
                     tok_per_commit, cost_per_rlog, tool_per_turn, role["compactions"], role["clears"]))
    body = ""
    for r, tokens, commits, ins, dele, rlog, tpc, cpr, tpt, compactions, clears in rows:
        body += (
            f"<tr><td class='role-cell' style='--role-color:{data['role_color'].get(r,'#888')}'>{_esc(r)}</td>"
            f"<td>{_fmt_int(tokens)}</td>"
            f"<td>{commits}</td>"
            f"<td>+{ins}/−{dele}</td>"
            f"<td>{rlog}</td>"
            f"<td>{_fmt_int(tpc) if tpc else '—'}</td>"
            f"<td>{_fmt_money(cpr) if cpr else '—'}</td>"
            f"<td>{tpt:.1f}</td>"
            f"<td>{compactions}</td>"
            f"<td>{clears}</td></tr>"
        )
    return (
        "<table class='data-table'><thead><tr>"
        "<th>Роль</th><th>Токены</th><th>Коммиты</th><th>+/− строк</th>"
        "<th>Записей в RESEARCH_LOG</th><th>Токенов/коммит</th><th>$/запись в логе</th>"
        "<th>Вызовов инстр./реплику</th><th>Compact-ов</th><th>/clear-ов</th>"
        "</tr></thead><tbody>" + body + "</tbody></table>"
    )


def table_chats(data: dict) -> str:
    body = ""
    chats = sorted(
        [c for c in data["chats"] if c.get("found")],
        key=lambda c: c.get("first_ts") or "",
    )
    for c in chats:
        body += (
            f"<tr><td class='role-cell' style='--role-color:{data['role_color'].get(c['role'],'#888')}'>{_esc(c['role'])}</td>"
            f"<td>{_esc(c['nick'])}</td>"
            f"<td><code>{_esc(c['chat_id'][:8])}…</code></td>"
            f"<td>{_esc((c.get('first_ts') or '')[:16].replace('T',' '))}</td>"
            f"<td>{_esc((c.get('last_ts') or '')[:16].replace('T',' '))}</td>"
            f"<td>{_fmt_h(c['duration_hours'])}</td>"
            f"<td>{c['user_turns']}</td>"
            f"<td>{c['tool_calls_total']}</td>"
            f"<td>{_fmt_int(c['tokens_grand_total'])}</td>"
            f"<td>{_fmt_money(c['cost_usd'])}</td>"
            f"<td>{c['compactions']}</td>"
            f"<td>{c.get('clears', 0)}</td></tr>"
        )
    missing = [c for c in data["chats"] if not c.get("found")]
    miss_note = ""
    if missing:
        names = ", ".join(f"{c['role']}/{c['nick']}" for c in missing)
        miss_note = f"<p class='callout'>Транскрипт не найден локально: {_esc(names)} " \
                    f"(вероятно, чат Cursor или другого клиента, не Claude Code CLI).</p>"
    return (
        "<table class='data-table'><thead><tr>"
        "<th>Роль</th><th>Ник</th><th>Chat ID</th><th>Начало (UTC)</th><th>Конец (UTC)</th>"
        "<th>Длительность</th><th>Реплик польз.</th><th>Вызовов инстр.</th><th>Токены</th>"
        "<th>~Стоимость</th><th>Compact-ов</th><th>/clear-ов</th>"
        "</tr></thead><tbody>" + body + "</tbody></table>" + miss_note
    )


def top_tools_table(data: dict) -> str:
    agg = {}
    for c in data["chats"]:
        for name, cnt in c.get("tool_calls", {}).items():
            agg[name] = agg.get(name, 0) + cnt
    top = sorted(agg.items(), key=lambda x: -x[1])[:12]
    body = "".join(f"<tr><td>{_esc(name)}</td><td>{cnt}</td></tr>" for name, cnt in top)
    return f"<table class='data-table'><thead><tr><th>Инструмент</th><th>Вызовов (все чаты)</th></tr></thead><tbody>{body}</tbody></table>"


# ---------------------------------------------------------------------------
# Главный шаблон
# ---------------------------------------------------------------------------

CSS = """
<style>
:root{
  --bg:#f7f4ee; --panel:#ffffff; --border:#e4ded2; --fg:#1c1b18; --muted:#6d6656;
  --accent:#b8650f; --accent-fg:#ffffff; --code-bg:#efe9db;
  --grid-line:#e4ded2; --callout-bg:#fdf1dc; --callout-fg:#8a5410; --callout-border:#eaceA0;
  --ok:#2f8f5b; --off:#8a8374;
}
@media (prefers-color-scheme: dark){
  :root{ --bg:#0f1113; --panel:#171a1e; --border:#2a2e33; --fg:#e9ece7; --muted:#8b9199;
    --accent:#ffb545; --accent-fg:#241300; --code-bg:#1d2126;
    --grid-line:#262a2f; --callout-bg:#241c0c; --callout-fg:#f0c46b; --callout-border:#4a3a17;
    --ok:#59c98a; --off:#7a7f87; }
}
:root[data-theme="dark"]{ --bg:#0f1113; --panel:#171a1e; --border:#2a2e33; --fg:#e9ece7; --muted:#8b9199;
  --accent:#ffb545; --accent-fg:#241300; --code-bg:#1d2126;
  --grid-line:#262a2f; --callout-bg:#241c0c; --callout-fg:#f0c46b; --callout-border:#4a3a17;
  --ok:#59c98a; --off:#7a7f87; }
:root[data-theme="light"]{ --bg:#f7f4ee; --panel:#ffffff; --border:#e4ded2; --fg:#1c1b18; --muted:#6d6656;
  --accent:#b8650f; --accent-fg:#ffffff; --code-bg:#efe9db;
  --grid-line:#e4ded2; --callout-bg:#fdf1dc; --callout-fg:#8a5410; --callout-border:#eaceA0;
  --ok:#2f8f5b; --off:#8a8374; }

*{box-sizing:border-box;}
body{background:var(--bg); color:var(--fg);}
.mono{font-family:ui-monospace,SFMono-Regular,"Cascadia Code",Consolas,"Liberation Mono",monospace;}
.tnum{font-variant-numeric:tabular-nums;}
.kpi-wrap{
  max-width:1020px; margin:0 auto; padding:28px 18px 72px;
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif; line-height:1.55;
}
h1{
  font-family:ui-monospace,SFMono-Regular,"Cascadia Code",Consolas,monospace;
  font-size:1.5rem; font-weight:700; margin:0 0 4px; letter-spacing:-.01em; text-wrap:balance;
}
h2{
  font-family:ui-monospace,SFMono-Regular,"Cascadia Code",Consolas,monospace;
  font-size:.78rem; font-weight:700; text-transform:uppercase; letter-spacing:.08em;
  color:var(--muted); margin:0 0 12px;
}
.panel + .panel, .panel-group + h2, h2 + .panel, h2 + .panel-group{margin-top:6px;}
section{margin-top:34px;}
.subtitle{color:var(--muted); font-size:.88rem; margin:0 0 4px;}
.subtitle.tight{margin-bottom:14px;}
.callout{
  background:var(--callout-bg); color:var(--callout-fg); border:1px solid var(--callout-border);
  border-radius:6px; padding:11px 15px; font-size:.86rem; margin:16px 0 0;
}
.callout code{background:transparent; padding:0; color:inherit; border-bottom:1px dotted currentColor;}

.cards{display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px;}
.card{
  background:var(--panel); border:1px solid var(--border); border-top:3px solid var(--role-color);
  border-radius:8px; padding:14px 14px 12px; display:flex; flex-direction:column; gap:2px;
}
.card-head{display:flex; align-items:center; justify-content:space-between; margin-bottom:6px;}
.card .role{font-weight:700; font-size:.92rem; font-family:ui-monospace,monospace;}
.chip{
  font-size:.62rem; font-weight:700; text-transform:uppercase; letter-spacing:.04em;
  padding:2px 7px; border-radius:99px; line-height:1.6;
}
.chip.on{background:color-mix(in srgb, var(--ok) 18%, transparent); color:var(--ok);}
.chip.off{background:color-mix(in srgb, var(--off) 18%, transparent); color:var(--off);}
.card .metric{font-size:1.32rem; font-weight:700; font-family:ui-monospace,monospace; margin-top:2px;}
.card .metric.small{font-size:1rem; margin-top:9px;}
.card .label{font-size:.68rem; color:var(--muted); text-transform:uppercase; letter-spacing:.04em;}

.panel{
  background:var(--panel); border:1px solid var(--border); border-radius:8px;
  padding:16px 18px 14px;
}
.panel + .panel{margin-top:12px;}

.chart-label{fill:var(--fg); font-size:11.5px; font-family:ui-monospace,monospace;}
.chart-legend{fill:var(--muted); font-size:10.5px; font-family:ui-monospace,monospace;}
.grid-line{stroke:var(--grid-line); stroke-width:1;}

.table-scroll{overflow-x:auto; border:1px solid var(--border); border-radius:8px;}
.data-table{width:100%; border-collapse:collapse; font-size:.82rem; min-width:640px;}
.data-table th{
  text-align:left; color:var(--muted); font-weight:600; padding:8px 10px;
  border-bottom:1px solid var(--border); font-size:.68rem; text-transform:uppercase; letter-spacing:.04em;
  white-space:nowrap; position:sticky; top:0; background:var(--panel);
}
.data-table td{padding:7px 10px; border-bottom:1px solid var(--border); font-variant-numeric:tabular-nums;}
.data-table tbody tr:hover{background:color-mix(in srgb, var(--accent) 6%, transparent);}
.role-cell{border-left:4px solid var(--role-color); font-weight:700; font-family:ui-monospace,monospace;}
code{background:var(--code-bg); padding:1px 6px; border-radius:4px; font-size:.85em; font-family:ui-monospace,monospace;}

.scenario{
  background:var(--panel); border:1px solid var(--border); border-left:3px solid var(--accent);
  border-radius:0 8px 8px 0; padding:11px 15px; display:flex; flex-direction:column; gap:3px;
}
.scenario b{color:var(--fg); font-size:.9rem;}
.scenario span{color:var(--muted); font-size:.85rem;}
.scenario-grid{display:grid; gap:8px;}

ul.method{padding-left:1.1em; display:flex; flex-direction:column; gap:8px; font-size:.88rem;}
ul.method li::marker{color:var(--accent);}

footer{color:var(--muted); font-size:.76rem; margin-top:44px; border-top:1px solid var(--border); padding-top:12px;}
a{color:var(--accent);}
:focus-visible{outline:2px solid var(--accent); outline-offset:2px;}
</style>
"""

SCENARIOS = [
    ("Где горит бюджет?", "Раздел «Токены по ролям» и «Стоимость по ролям» — сразу видно, какая роль "
     "тяжелее всех по токенам/деньгам. Если это не совпадает с ощущением важности роли — повод "
     "пересмотреть промт или разбить задачу на более мелкие шаги."),
    ("Кто продуктивнее на единицу токенов?", "Таблица «Эффективность»: токенов на коммит и $ на "
     "запись в RESEARCH_LOG. Низкое число = дёшево произведённый результат; высокое = либо сложная "
     "задача, либо чат «буксует» (много рассуждений, мало фиксируемого результата)."),
    ("Не залипает ли чат в лимит контекста?", "Столбец «Compact-ов» в детальной таблице — сколько раз "
     "чат упирался в `/compact`. Много compact-ов у одной роли = её промт/файлы координации слишком "
     "тяжёлые для повторного чтения на старте, стоит сжать CHARTER/BOARD или чаще делать хэндоффы."),
    ("Были ли параллельные (конфликтующие) чаты одной роли?", "График «Хронология жизни чатов» — "
     "если два бара одной роли перекрываются по времени, это ровно та ситуация, от которой защищает "
     "CHARTER §12 (зафиксированный прецедент - L1/L2, 29-30 июля)."),
    ("Сравнить два чата одной роли (напр. A1 vs A2)", "Детальная таблица отсортирована по времени "
     "начала — соседние строки одной роли сравнимы напрямую: сколько токенов/времени/коммитов дал "
     "каждый конкретный чат-инстанс."),
    ("Какие инструменты самые дорогие по вызовам?", "Таблица «Топ инструментов» — если Bash/Read "
     "доминируют с большим отрывом, возможно, часть работы стоит закрыть более точечными "
     "инструментами или скриптами вместо интерактивного перебора."),
]


def render(data: dict) -> str:
    role_order = data["role_order"]
    role_color = data["role_color"]
    roles = data["roles"]

    cards = ""
    for r in role_order:
        role = roles[r]
        chip = "<span class='chip on'>active</span>" if role["is_active"] else "<span class='chip off'>stopped</span>"
        cards += (
            f"<div class='card' style='--role-color:{role_color.get(r,'#888')}'>"
            f"<div class='card-head'><span class='role'>{_esc(r)}</span>{chip}</div>"
            f"<div class='metric tnum'>{_fmt_int(role['tokens_grand_total'])}</div>"
            f"<div class='label'>токенов &middot; {role['chats']} чат(ов)</div>"
            f"<div class='metric small tnum'>{_fmt_money(role['cost_usd'])}</div>"
            f"<div class='label'>оценка стоимости</div>"
            f"<div class='metric small tnum'>{_fmt_h(role['duration_hours'])}</div>"
            f"<div class='label'>суммарное время жизни чатов</div>"
            "</div>"
        )

    tokens_chart = svg_stacked_tokens(roles, role_order, role_color)
    cost_items = [(r, roles[r]["cost_usd"]) for r in role_order if roles[r]["cost_usd"] > 0]
    cost_items.sort(key=lambda x: -x[1])
    cost_chart = svg_simple_bars(cost_items, role_color, fmt=_fmt_money)
    commits_items = [(r, data["git"][r]["commits"]) for r in role_order if data["git"][r]["commits"] > 0]
    commits_items.sort(key=lambda x: -x[1])
    commits_chart = svg_simple_bars(commits_items, role_color, fmt=lambda v: str(int(v)))
    gantt_chart = svg_gantt(data["chats"], role_color)

    scenarios_html = "".join(
        f"<div class='scenario'><b>{_esc(q)}</b><span>{_esc(a)}</span></div>" for q, a in SCENARIOS
    )

    generated = data["generated_at"][:19].replace("T", " ") + " UTC"

    return f"""
<title>KPI-отчёт: THz-Unified-Optimizer</title>
{CSS}
<div class="kpi-wrap">
  <h1>Приборная панель координации &middot; THz-Unified-Optimizer</h1>
  <div class="subtitle">сгенерировано {_esc(generated)} &middot;
    источник — <code>NICKS.md</code> + локальные транскрипты Claude Code + <code>git log</code> +
    <code>RESEARCH_LOG.md</code></div>

  <div class="callout">
    ⚠ Оценка стоимости ($) — ОРИЕНТИРОВОЧНАЯ (тарифы для модельной линейки этого проекта официально не
    подтверждены, перенесены из ближайшего известного прайс-листа). Годится для СРАВНЕНИЯ ролей между
    собой, не для сверки со счётом. Токены — точные (взяты из полей <code>usage</code> транскриптов).
  </div>

  <section>
    <h2>Сводка по ролям</h2>
    <div class="cards">{cards}</div>
  </section>

  <section>
    <h2>Токены по ролям — input / output / cache-write / cache-read</h2>
    <div class="panel">{tokens_chart}</div>
  </section>

  <section>
    <h2>Оценка стоимости по ролям</h2>
    <div class="panel">{cost_chart}</div>
  </section>

  <section>
    <h2>Коммиты по ролям — git log, тег [РОЛЬ]</h2>
    <div class="panel">{commits_chart}</div>
  </section>

  <section>
    <h2>Хронология жизни чатов</h2>
    <p class="subtitle tight">Каждый бар — один физический чат (ник) от первого до последнего
      сообщения. Перекрывающиеся бары одной роли = коллизия (см. CHARTER §12, прецедент L1/L2).</p>
    <div class="panel">{gantt_chart}</div>
  </section>

  <section>
    <h2>Эффективность по ролям</h2>
    <div class="table-scroll">{table_efficiency(data)}</div>
  </section>

  <section>
    <h2>Топ инструментов по вызовам — все чаты</h2>
    <div class="table-scroll">{top_tools_table(data)}</div>
  </section>

  <section>
    <h2>Детализация по чатам</h2>
    <div class="table-scroll">{table_chats(data)}</div>
  </section>

  <section>
    <h2>Как читать этот отчёт — сценарии использования</h2>
    <div class="scenario-grid">{scenarios_html}</div>
  </section>

  <section>
    <h2>Методология и оговорки</h2>
  <ul class="method">
    <li><b>Токены</b> — сумма полей <code>usage.input_tokens/output_tokens/cache_creation_input_tokens/
      cache_read_input_tokens</code> по всем assistant-сообщениям транскрипта (модель <code>&lt;synthetic&gt;</code>
      исключена — это служебные, не биллингуемые записи).</li>
    <li><b>«Реплик пользователя»</b> — user-записи, содержащие хотя бы один текстовый блок; записи,
      целиком состоящие из <code>tool_result</code> (автоматическая обратная связь инструментов), не
      считаются репликой человека.</li>
    <li><b>Compact-ы</b> — вхождения <code>&lt;command-name&gt;compact&lt;/command-name&gt;</code> в
      тексте реплик — то есть сколько раз ЭТОТ чат достиг предела контекста и был сжат.</li>
    <li><b>Коммиты/строки</b> — <code>git log</code>, коммиты с subject, начинающимся на
      <code>[РОЛЬ]</code>; для B в этом репозитории таких коммитов 0 (её первая работа попала в общий
      бутстрап-коммит без тега — известный факт, см. <code>sessions/B.md</code>), это не значит «ноль
      работы».</li>
    <li><b>RESEARCH_LOG</b> — считает только записи в общем файле; роли, ведущие числа в других файлах
      (D → <code>DRAFT.md</code>, L → <code>litrev/answers/**</code>, B → <code>CHANGELOG_model_core.md</code>)
      здесь недоучтены — раздел показывает нижнюю границу, не полный объём.</li>
    <li><b>Sidechain/субагенты</b> — если чат использовал инструмент Agent (подзадачи), их токены могут
      быть в ОТДЕЛЬНЫХ файлах транскриптов, не привязанных напрямую к нику в <code>NICKS.md</code> —
      скрипт их не подтягивает (нет надёжного способа сопоставить subagent-транскрипт с родительским
      чатом без дополнительных данных).</li>
  </ul>
  </section>

  <footer>
    Сгенерировано <code>coordination/tools/kpi_stats.py</code>. Перезапуск:
    <code>python coordination/tools/kpi_stats.py</code> из корня репозитория. Зона ORCH (CHARTER §3).
  </footer>
</div>
"""
