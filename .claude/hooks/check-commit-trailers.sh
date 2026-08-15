#!/usr/bin/env bash
# Проверка разметки коммита (CHARTER §4).
#
# Зачем: за неделю 09–15.08 трейлер `Session:` оказался нераспознаваемым git
# в 10 коммитах из 29 — текст был на месте, но git его не видел. Две причины,
# обе механические: пустая строка между трейлерами и подписью, либо `Reason:`,
# перенесённый на вторую строку без отступа (тогда весь блок перестаёт быть
# блоком трейлеров). Правило в уставе от этого не спасало — проверка спасает.
#
# Как: не разбираем текст команды (heredoc надёжно не распарсить), а спрашиваем
# у самого git, видит ли ОН трейлеры в получившемся коммите. Парсер git —
# единственный авторитет в этом вопросе.
#
# Событие: PostToolUse / Bash. Вердикт `block` возвращает причину модели,
# ход продолжается — коммит чинится через `git commit --amend`.
set -uo pipefail

input=$(cat)

# Интересует только git commit. Ищем подстрокой: в JSON она не экранируется.
case "$input" in
  *"git commit"*) ;;
  *) exit 0 ;;
esac

git rev-parse --git-dir >/dev/null 2>&1 || exit 0

# Если коммит не состоялся (команда упала), HEAD остался старым — не трогаем.
now=$(date +%s)
made=$(git log -1 --format=%ct 2>/dev/null || echo 0)
[ $((now - made)) -gt 120 ] && exit 0

session=$(git log -1 --format='%(trailers:key=Session,valueonly)' 2>/dev/null | tr -d '[:space:]')
reason=$(git log -1 --format='%(trailers:key=Reason,valueonly)'  2>/dev/null | tr -d '[:space:]')

miss=""
[ -z "$session" ] && miss="Session:"
[ -z "$reason" ]  && miss="${miss:+$miss и }Reason:"
[ -z "$miss" ] && exit 0

subject=$(git log -1 --format=%s 2>/dev/null | sed 's/[\\"]/ /g' | cut -c1-70)

printf '{"decision":"block","reason":"Коммит «%s»: git НЕ РАСПОЗНАЁТ трейлер %s — текст может быть на месте, но парсер его не видит. Причина почти всегда одна из двух: (1) пустая строка между трейлерами и подписью Co-Authored-By, (2) Reason: перенесён на вторую строку. Трейлеры обязаны быть ПОСЛЕДНИМ абзацем без пустых строк внутри, Reason: — в одну строку (CHARTER §4). Почини: git commit --amend -F - с корректным сообщением, затем проверь: git log -1 --pretty=%%(trailers:key=Session,valueonly)"}\n' \
  "$subject" "$miss"
exit 0
