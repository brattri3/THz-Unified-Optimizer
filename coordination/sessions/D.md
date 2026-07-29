# Сессия D — метка намерений / замок

> Обновляй при СТАРТЕ, смене задачи, «пульсе» и остановке. Перед правкой ОБЩЕГО файла — проверь
> чужие sessions/*.md на пересечение. Устаревший heartbeat (>~1ч) = сессия считается неактивной.

- **Статус:** running (запущена владельцем 2026-07-29 13:19)
- **Специализация:** статья + **образовательный слой** (CHARTER §3, §11 / OWNERSHIP.md)
- **Текущая задача:** D1 (обзор для статьи) — **done**; D2 (запрос к L) — **оформлен**;
  далее D3 (DRAFT) и D4 (primers)
- **План (очередь задач D):**
  - **D0** регистрация + `research/paper/STATE.md` + принятие seed `problem_statement.tex` — done
  - **D1** обзор литературы ДЛЯ СТАТЬИ из имеющихся источников → `research/paper/RELATED_WORK.md` — done
  - **D2** структурированный запрос `ОТ D К L` в `HANDOFFS.md` по пробелам D1 — done
  - **D3** `research/paper/DRAFT.md` — Abstract / Methods / Results / Novelty / Limitations
  - **D4** `research/paper/primers/` — методички: матрицы Джонса и поляризация; модель Бланко и
    эквивалентные схемы; модель Друде; рассеяние на неоднородностях (ν^γ); нелинейный МНК и
    Левенберг–Марквардт; правдоподобие и AIC; информация Фишера / идентифицируемость;
    бутстрап и профильное правдоподобие
  - **D5** `GLOSSARY.md` — ведётся сквозным образом (v0.1 создан вместе с D1)
- **Файлы, что трогаю (замок):**
  `research/paper/**` (STATE.md, RELATED_WORK.md, GLOSSARY.md, DRAFT.md, primers/**) —
  эксклюзивная зона D; плюс СВОИ строки/записи в общих `coordination/BOARD.md`,
  `coordination/ACTIVITY.md`, `coordination/HANDOFFS.md`, этот файл.
- **Заявлено до (время):** до остановки сессии владельцем; heartbeat обновляется каждой итерацией /loop
- **Heartbeat:** 2026-07-29 13:19

## Границы (соблюдаю)
- Первичный литературный поиск НЕ веду — заказываю у **L** через `HANDOFFS.md`
  (OWNERSHIP «Граница D↔L»: L = поиск/извлечение и обзор ПОЛЯ, D = синтез и обзор ДЛЯ СТАТЬИ).
- Новые расчёты/фиты НЕ делаю — запрос к **A** (физика/гипотезы) или **B** (численное ядро).
- Пишу только в `research/paper/**`; `research/papers/**`, `research/results/**`,
  `research/two_wgp/**`, `research/SYNTHESIS.md`, `research/logs/RESEARCH_LOG.md` — только чтение.

## История намерений (append)
- [2026-07-27] seed-инициализация метки.
- [2026-07-29, ORCH] Уточнён профиль D: добавлен образовательный мандат (владелец — не специалист в
  предметной математике/физике, тексты D обязаны обучать; CHARTER §11). Промт запуска обновлён в
  LAUNCH_PROMPTS.md. В зону D положен seed `problem_statement.tex` (постановка задачи). D пока не
  запущена — стартует владельцем готовой строкой.
- [2026-07-29 13:19, D] **Старт.** Прочитаны CHARTER (вкл. §11), ORCH_BRIEF, OWNERSHIP, BOARD,
  HANDOFFS, SYNTHESIS, seed `problem_statement.tex`, а также источники L: `research/papers/INDEX.md`,
  `litrev/RETROSPECTIVE_REVIEW.md`, `litrev/LITERATURE_REVIEW.md` (перечень corpus/*.json и txt/*.txt).
  Заведены `research/paper/STATE.md`, `RELATED_WORK.md` (D1), `GLOSSARY.md`. Оформлен запрос
  `ОТ D К L` (D2). Конфликтов зон нет: A работает в `research/two_wgp/**` (A0/A6), L — в
  `research/papers/**`; я пишу только в `research/paper/**`.
