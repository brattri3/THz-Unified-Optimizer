# BOARD — живой статус сессий

Каждая сессия обновляет СВОЮ строку: при старте, смене задачи, «пульсе» (раз в несколько шагов) и
остановке. Формат времени — ISO (локальное). Статус: `running` / `idle` / `blocked` / `stopped`.

| Сессия | Статус | Текущая задача | План (кратко) | Обновлено (heartbeat) |
|---|---|---|---|---|
| **A** модели/идентифицируемость/гипотезы + софт затухания | idle/stale (замок снят ORCH 17:16) | **A6_angle_symmetry — В РАБОТЕ, НЕ ЗАВЕРШЕНА** | Скрипты написаны (`a6_angle_symmetry.py` 452стр, `exp_builder.py` 192стр, незакоммичены), но НЕТ результата: нет json в results/two_wgp, нет записи в RESEARCH_LOG, state.json=pending. Сессия встала ~08:24. Возобновить: запустить A6, залогировать числа, закрыть задачу. Далее → A4 (API для C) → A3 → A5 | 2026-07-28 08:24 |
| **B** ядро/рефактор | idle | model_core доведён (итер. 1–2), бит-в-бит | далее по запросу: two_wgp на ядро (нужен хэндофф A), перегенерация устаревших артефактов | 2026-07-28 |
| **C** app | idle | **C6_v02_acceptance** — UI переведён на EN (прибор в Китае, решение владельца, ACTIVITY 21:55), физика бит-в-бит; Win7-комплект + чек-лист ред.2 распечатан | ждём прогон владельца на Win7-ПК у спектрометра + обратную связь → флаг одобрения v0.2; далее C2 (ждёт A4), C3–C5. Ветка `session/C` (`../TUO-C`, HEAD `713f807`) готова к мержу ORCH | 2026-07-28 21:55 |
| **D** статья + образоват. слой | seed (не запускалась; промт готов) | — | При старте: принять seed `problem_statement.tex`, DRAFT + `primers/` (методички) + GLOSSARY; педагог. мандат (CHARTER §11). Расчёты — через HANDOFFS к A/B | 2026-07-27 (seed) |
| **L** литература/обзор + метапоиск | running | **L1_review_maintain** | Поддержка и углубление ретро-обзора (`litrev/RETROSPECTIVE_REVIEW.md`), граф цитирований через `metasearch.py` (OpenAlex/Crossref/Unpaywall/arXiv), ГОСТ-библиография с DOI, снабжение D (статья) и A (гипотезы) источниками. Отслеживание конкурентов (Karimi 2026, PureWave). Зона: `research/papers/**` (litrev свёрнут в наш репо ORCH 2026-07-28). Запущена с санкции владельца (CHARTER §8) | 2026-07-28 |
| **ORCH** control | running (проход 2026-07-28 17:16) | свод + git-гигиена | Снят устаревший замок A (heartbeat 08:24). Конфликтов зон нет. Открытые HANDOFFS: B→A blanco_t (open), C→A A4-crosscheck (open, ждёт A4), B→A build_experiment (taken). Флаг владельцу: незакоммиченная работа A/C/ORCH на main | 2026-07-28 17:16 |

## Легенда владения (полностью — OWNERSHIP.md)
A→`research/two_wgp/`+`research/hypotheses/` · B→`research/experiments/model_core.py`+shared-code · C→`attenuator_app/` ·
D→`research/paper/` · L→`research/papers/` (вкл. свёрнутый `research/papers/litrev/`) · ORCH→`coordination/`,`CLAUDE.md`,main `research/state.json`.

## Активные «замки» намерений
См. `coordination/sessions/*.md`. Перед правкой общего файла — проверь их.
