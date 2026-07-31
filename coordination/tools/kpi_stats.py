"""
KPI-статистика мультисессионной работы над THz-Unified-Optimizer (зона ORCH, CHARTER §12/NICKS.md).

Что делает:
1. Читает реестр чатов `coordination/NICKS.md` (ник -> роль -> chat id).
2. Для каждого chat id ищет транскрипт Claude Code в `~/.claude/projects/**/<id>.jsonl`
   и считает по нему токены, вызовы инструментов, длительность, число compact-событий.
3. Кросс-сверяет с `git log` (коммиты/строки по тегу `[<РОЛЬ>]`) и с
   `research/logs/RESEARCH_LOG.md` (записи по тегу `[<РОЛЬ>]`).
4. Пишет сырые данные в `coordination/reports/kpi_data.json` и рендерит
   самодостаточный HTML-отчёт `coordination/reports/kpi_report.html`.

Запуск (из корня репозитория или откуда угодно - пути абсолютные):
    python coordination/tools/kpi_stats.py

Ничего не разбирается о РЕЗУЛЬТАТАХ исследования - это чисто метрики РАБОТЫ чатов
(токены, время, коммиты). Наука/выводы А/Б/В экспериментов сюда не тянутся.

⚠ ОЦЕНКА СТОИМОСТИ ($) - ОРИЕНТИРОВОЧНАЯ. Таблица тарифов ниже перенесена с последнего
известного прайс-листа Anthropic API (не подтверждена для модельной линейки "5" из этого
проекта) и нужна ТОЛЬКО для относительного сравнения сессий между собой, а не для сверки
с реальным счётом (если аккаунт вообще на pay-per-token, а не на подписке). Обновите
таблицу PRICING_USD_PER_MTOK, если появится точный прайс.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
COORD = REPO / "coordination"
NICKS_MD = COORD / "NICKS.md"
RESEARCH_LOG = REPO / "research" / "logs" / "RESEARCH_LOG.md"
REPORTS_DIR = COORD / "reports"
CLAUDE_PROJECTS = Path.home() / ".claude" / "projects"

ROLE_ORDER = ["ORCH", "A", "B", "C", "D", "L"]
ROLE_COLOR = {
    "ORCH": "#8890a0",
    "A": "#4fb3e0",
    "B": "#2fbf98",
    "C": "#e2933f",
    "D": "#b876e0",
    "L": "#e0645f",
}

# $ за 1M токенов. ОРИЕНТИРОВОЧНО - см. предупреждение в докстринге модуля.
PRICING_USD_PER_MTOK = {
    "claude-opus-5":    {"input": 15.0, "output": 75.0, "cache_write": 18.75, "cache_read": 1.50},
    "claude-sonnet-5":  {"input": 3.0,  "output": 15.0, "cache_write": 3.75,  "cache_read": 0.30},
    "claude-fable-5":   {"input": 3.0,  "output": 15.0, "cache_write": 3.75,  "cache_read": 0.30},
    "claude-haiku-4-5": {"input": 0.8,  "output": 4.0,  "cache_write": 1.0,  "cache_read": 0.08},
}
DEFAULT_RATE = PRICING_USD_PER_MTOK["claude-sonnet-5"]


def parse_nicks(path: Path) -> list[dict]:
    rows = []
    id_re = re.compile(r"`([0-9a-fA-F-]{20,36})`")
    nick_re = re.compile(r"\*\*([A-Za-z0-9]+)\*\*")
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line.startswith("|") or line.startswith("|---"):
            continue
        cols = [c.strip() for c in line.split("|")][1:-1]
        if len(cols) < 6:
            continue
        nick_m = nick_re.search(cols[0])
        id_m = id_re.search(cols[2])
        if not (nick_m and id_m):
            continue
        rows.append({
            "nick": nick_m.group(1),
            "role": cols[1],
            "chat_id": id_m.group(1),
            "status_note": cols[5],
            "is_active": "active" in cols[5].lower(),
        })
    return rows


def find_transcript(chat_id: str) -> Path | None:
    matches = list(CLAUDE_PROJECTS.glob(f"*/{chat_id}.jsonl"))
    return matches[0] if matches else None


def parse_ts(ts: str) -> datetime | None:
    if not ts:
        return None
    try:
        return datetime.fromisoformat(ts.replace("Z", "+00:00"))
    except ValueError:
        return None


@dataclass
class ChatStats:
    nick: str
    role: str
    chat_id: str
    found: bool = False
    project_dir: str = ""
    first_ts: str | None = None
    last_ts: str | None = None
    duration_hours: float = 0.0
    user_entries: int = 0
    user_turns: int = 0          # user-записи, содержащие хоть один текстовый блок (не чистый tool_result)
    assistant_entries: int = 0
    compactions: int = 0
    clears: int = 0
    tool_calls: Counter = field(default_factory=Counter)
    effort_counts: Counter = field(default_factory=Counter)
    tokens_by_model: dict = field(default_factory=dict)  # model -> {input,output,cache_write,cache_read}

    def tokens_total(self) -> dict:
        out = {"input": 0, "output": 0, "cache_write": 0, "cache_read": 0}
        for m in self.tokens_by_model.values():
            for k in out:
                out[k] += m.get(k, 0)
        return out

    def cost_usd(self) -> float:
        total = 0.0
        for model, toks in self.tokens_by_model.items():
            rate = PRICING_USD_PER_MTOK.get(model, DEFAULT_RATE)
            total += toks.get("input", 0) / 1e6 * rate["input"]
            total += toks.get("output", 0) / 1e6 * rate["output"]
            total += toks.get("cache_write", 0) / 1e6 * rate["cache_write"]
            total += toks.get("cache_read", 0) / 1e6 * rate["cache_read"]
        return total


def analyze_transcript(path: Path) -> ChatStats | None:
    stats = ChatStats(nick="", role="", chat_id=path.stem, found=True, project_dir=path.parent.name)
    first_dt = last_dt = None
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            ts = parse_ts(obj.get("timestamp"))
            if ts:
                if first_dt is None or ts < first_dt:
                    first_dt = ts
                if last_dt is None or ts > last_dt:
                    last_dt = ts
            t = obj.get("type")
            if t == "user":
                stats.user_entries += 1
                msg = obj.get("message", {})
                content = msg.get("content")
                has_text = False
                text_acc = ""
                if isinstance(content, str):
                    has_text = bool(content.strip())
                    text_acc = content
                elif isinstance(content, list):
                    for c in content:
                        if isinstance(c, dict) and c.get("type") == "text" and c.get("text", "").strip():
                            has_text = True
                            text_acc += c.get("text", "")
                if has_text:
                    stats.user_turns += 1
                    if "<command-name>/compact</command-name>" in text_acc or \
                       "<command-message>compact</command-message>" in text_acc:
                        stats.compactions += 1
                    if "<command-name>/clear</command-name>" in text_acc or \
                       "<command-message>clear</command-message>" in text_acc:
                        stats.clears += 1
            elif t == "assistant":
                stats.assistant_entries += 1
                msg = obj.get("message", {})
                model = msg.get("model", "unknown")
                effort = obj.get("effort")
                if effort:
                    stats.effort_counts[effort] += 1
                usage = msg.get("usage")
                if usage and model and model != "<synthetic>":
                    bucket = stats.tokens_by_model.setdefault(
                        model, {"input": 0, "output": 0, "cache_write": 0, "cache_read": 0}
                    )
                    bucket["input"] += usage.get("input_tokens", 0)
                    bucket["output"] += usage.get("output_tokens", 0)
                    bucket["cache_write"] += usage.get("cache_creation_input_tokens", 0)
                    bucket["cache_read"] += usage.get("cache_read_input_tokens", 0)
                content = msg.get("content")
                if isinstance(content, list):
                    for c in content:
                        if isinstance(c, dict) and c.get("type") == "tool_use":
                            stats.tool_calls[c.get("name", "?")] += 1
    stats.first_ts = first_dt.isoformat() if first_dt else None
    stats.last_ts = last_dt.isoformat() if last_dt else None
    if first_dt and last_dt:
        stats.duration_hours = round((last_dt - first_dt).total_seconds() / 3600, 2)
    return stats


def git_log_subjects() -> list[tuple[str, str]]:
    out = subprocess.run(
        ["git", "log", "--pretty=format:%H\x01%s"],
        cwd=REPO, capture_output=True, text=True, check=True,
        encoding="utf-8", errors="replace",
    ).stdout
    rows = []
    for line in out.splitlines():
        if "\x01" in line:
            h, s = line.split("\x01", 1)
            rows.append((h, s))
    return rows


def git_stats_by_tag() -> dict[str, dict]:
    subjects = git_log_subjects()
    by_tag: dict[str, dict] = {r: {"commits": 0, "insertions": 0, "deletions": 0} for r in ROLE_ORDER}
    for h, subj in subjects:
        m = re.match(r"^\[([A-Za-z]+)\]", subj)
        if not m:
            continue
        tag = m.group(1)
        if tag not in by_tag:
            continue
        shortstat = subprocess.run(
            ["git", "show", "--shortstat", "--format=", h],
            cwd=REPO, capture_output=True, text=True, check=True,
            encoding="utf-8", errors="replace",
        ).stdout.strip()
        ins = dele = 0
        mi = re.search(r"(\d+) insertion", shortstat)
        md = re.search(r"(\d+) deletion", shortstat)
        if mi:
            ins = int(mi.group(1))
        if md:
            dele = int(md.group(1))
        by_tag[tag]["commits"] += 1
        by_tag[tag]["insertions"] += ins
        by_tag[tag]["deletions"] += dele
    return by_tag


def research_log_counts() -> dict[str, int]:
    if not RESEARCH_LOG.exists():
        return {r: 0 for r in ROLE_ORDER}
    text = RESEARCH_LOG.read_text(encoding="utf-8")
    return {r: len(re.findall(rf"\[{re.escape(r)}\]", text)) for r in ROLE_ORDER}


def build_data() -> dict:
    nick_rows = parse_nicks(NICKS_MD)
    chats = []
    for row in nick_rows:
        nick, role, chat_id = row["nick"], row["role"], row["chat_id"]
        path = find_transcript(chat_id)
        if path is None:
            chats.append({
                "nick": nick, "role": role, "chat_id": chat_id, "found": False,
                "is_active": row["is_active"], "status_note": row["status_note"],
            })
            continue
        st = analyze_transcript(path)
        st.nick, st.role, st.chat_id = nick, role, chat_id
        tokens = st.tokens_total()
        chats.append({
            "nick": nick, "role": role, "chat_id": chat_id, "found": True,
            "is_active": row["is_active"], "status_note": row["status_note"],
            "project_dir": st.project_dir,
            "first_ts": st.first_ts, "last_ts": st.last_ts,
            "duration_hours": st.duration_hours,
            "user_turns": st.user_turns, "user_entries": st.user_entries,
            "assistant_entries": st.assistant_entries,
            "compactions": st.compactions,
            "clears": st.clears,
            "tool_calls": dict(st.tool_calls),
            "tool_calls_total": sum(st.tool_calls.values()),
            "effort_counts": dict(st.effort_counts),
            "tokens_by_model": st.tokens_by_model,
            "tokens_total": tokens,
            "tokens_grand_total": sum(tokens.values()),
            "cost_usd": round(st.cost_usd(), 4),
        })

    roles = {r: {
        "chats": 0, "duration_hours": 0.0, "user_turns": 0, "tool_calls_total": 0,
        "tokens_total": {"input": 0, "output": 0, "cache_write": 0, "cache_read": 0},
        "tokens_grand_total": 0, "cost_usd": 0.0, "compactions": 0, "clears": 0,
        "is_active": False,
    } for r in ROLE_ORDER}
    for c in chats:
        if c.get("is_active"):
            roles[c["role"]]["is_active"] = True
        if not c.get("found"):
            continue
        r = roles[c["role"]]
        r["chats"] += 1
        r["duration_hours"] += c["duration_hours"]
        r["user_turns"] += c["user_turns"]
        r["tool_calls_total"] += c["tool_calls_total"]
        for k in r["tokens_total"]:
            r["tokens_total"][k] += c["tokens_total"][k]
        r["tokens_grand_total"] += c["tokens_grand_total"]
        r["cost_usd"] += c["cost_usd"]
        r["compactions"] += c["compactions"]
        r["clears"] += c["clears"]
    for r in roles.values():
        r["cost_usd"] = round(r["cost_usd"], 4)
        r["duration_hours"] = round(r["duration_hours"], 2)

    git_stats = git_stats_by_tag()
    rlog = research_log_counts()

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "chats": chats,
        "roles": roles,
        "git": git_stats,
        "research_log": rlog,
        "pricing_table": PRICING_USD_PER_MTOK,
        "role_order": ROLE_ORDER,
        "role_color": ROLE_COLOR,
    }


def main():
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    data = build_data()
    json_path = REPORTS_DIR / "kpi_data.json"
    json_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Written: {json_path}")

    # рендер HTML - в отдельном модуле, чтобы держать логику сбора данных отдельно от вёрстки
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import kpi_report  # noqa: E402

    html_path = REPORTS_DIR / "kpi_report.html"
    html_path.write_text(kpi_report.render(data), encoding="utf-8")
    print(f"Written: {html_path}")


if __name__ == "__main__":
    main()
