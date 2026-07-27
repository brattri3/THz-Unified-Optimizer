# OWNERSHIP — карта владения путями

«Владеет» = имеет право ПИСАТЬ. Остальные — read-only; правка → запрос через `HANDOFFS.md` + `ACTIVITY.md`.

## Эксклюзивные зоны
| Путь / файл | Владелец | Прочим |
|---|---|---|
| `research/two_wgp/**` | A | read |
| `research/results/two_wgp/**` | A | read |
| `research/experiments/model_core.py`, `verify_model_core.py`, `refactor_ab_check.py` | B | read |
| `research/experiments/*.py` (общие фит-скрипты t*, model_m5, fit_lib) | **B** (как рефактор-владелец) | read; правки — через хэндофф B |
| `attenuator_app/**`, `docs/attenuator_app/**` | C | read |
| `research/paper/**` | D | read |
| `coordination/**`, `CLAUDE.md`, `research/state.json` (main), `research/RESEARCH_PLAN.md` | ORCH | read; предложения — в ACTIVITY |

## Общие (shared) файлы и правила доступа
| Файл | Правило |
|---|---|
| `research/logs/RESEARCH_LOG.md` | **append-only**, все сессии дописывают с тегом `[<ID>]`; не редактировать чужие записи |
| `research/hypotheses/HYPOTHESES.md` | пишет A (гипотезы) + ORCH; менять только по своей гипотезе, анонс в ACTIVITY |
| `data_pool/**` | **read-only** для всех агентов (гардрейл §6); данные добавляет только ВЛАДЕЛЕЦ-человек |
| `unified_optimizer/**` | **read-only** (продакшн; правка ломает базовую линию — спросить владельца) |
| `coordination/BOARD.md`, `sessions/<ID>.md`, `ACTIVITY.md` | каждая сессия пишет СВОЮ строку/файл/добавляет запись |

## Per-session state (не общий state.json!)
Каждая сессия ведёт СВОЙ state:
- A → `research/two_wgp/state.json`
- B → `research/experiments/CHANGELOG_model_core.md` (+ свой state при необходимости)
- C → `attenuator_app/STATE.md` (создать)
- D → `research/paper/STATE.md` (создать)
- Главная линия исследования → `research/state.json` (владелец ORCH / основной анализ).

## Конфликтные горячие точки (беречь особо)
- `fit_lib.py`, `model_core.py` — правит только B. A/C/D читают.
- `RESEARCH_LOG.md` — только append, тег сессии.
- `research/state.json` (main) vs `research/two_wgp/state.json` — не путать; каждая сессия свой.
- git `main` — при общем каталоге коммитят по очереди с анонсом; при worktree — мерж через ORCH.
