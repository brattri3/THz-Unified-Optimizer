#!/usr/bin/env bash
# Хук бюджета контекста холодного старта (PARADIGM_REVIEW_2026-08-23 §7 п.1).
#
# Зачем: калибровка "цена в CLAUDE.md" протухала молча (Д-2) — файлы росли,
# таблица нет, роль решала "читать/не читать" по неверной цифре. Правило
# "ужать до 1200 токенов" из COLDSTART_REDESIGN_2026-08-15 не проверялось
# нигде. Правила, встроенные в действие, держатся сами (CHARTER §4) — этот
# хук переводит "не забыть" в "проверяется", как уже сделано для трейлеров.
#
# Единица — байты UTF-8 (`wc -c`), не Unicode code points и не токены: это
# то, чем на самом деле измерена вся история проекта (roles/A.md — 9739
# байт ≈ 4870 ток. × 2.0 симв/ток из калибровки — сходится с байтами, не
# с числом символов). См. budget.json.
#
# Событие: SessionStart, только source=startup|clear — resume/compact файлы
# старта заново не читают, предупреждать там — чистый шум.
set -uo pipefail
input=$(cat)

case "$input" in
  *'"source":"startup"'*|*'"source": "startup"'*) ;;
  *'"source":"clear"'*|*'"source": "clear"'*) ;;
  *) exit 0 ;;
esac

root="${CLAUDE_PROJECT_DIR:-}"
[ -d "$root" ] || root=$(git rev-parse --show-toplevel 2>/dev/null)
[ -d "$root" ] || exit 0

budget="$root/.claude/hooks/budget.json"
[ -f "$budget" ] || exit 0

limit=$(grep -o '"limit_bytes"[[:space:]]*:[[:space:]]*[0-9]*' "$budget" | head -1 | grep -o '[0-9]*$')
[ -n "${limit:-}" ] || exit 0

over=""
for f in "$root"/coordination/roles/*.md; do
  [ -f "$f" ] || continue
  size=$(wc -c < "$f" | tr -d '[:space:]')
  [ "$size" -gt "$limit" ] || continue
  over="${over}${over:+, }$(basename "$f") ($size/$limit Б)"
done
[ -n "$over" ] || exit 0

printf '{"systemMessage":"Бюджет холодного старта превышен: %s. Лимит coordination/roles/*.md — %s байт (budget.json)."}\n' "$over" "$limit"
exit 0
