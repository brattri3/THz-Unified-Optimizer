# OWNERSHIP — карта владения путями

«Владеет» = имеет право ПИСАТЬ. Остальные — read-only; правка → запрос через `HANDOFFS.md` + `ACTIVITY.md`.

## Эксклюзивные зоны
| Путь / файл | Владелец | Прочим |
|---|---|---|
| `research/two_wgp/**` | A | read |
| `research/results/two_wgp/**` | A | read |
| `research/hypotheses/**` | A | read (см. санкцию владельца 2026-07-28, CHARTER §3) |
| **Численное ядро** (где бы ни лежало): `research/experiments/model_core.py`, `fit_lib.py`, `verify_model_core.py`, `refactor_ab_check.py` | **B** + роль ревьюера численной эквивалентности для правок ядра ЛЮБОЙ сессией | read |
| **Скрипты-эксперименты/гипотезы** `research/experiments/tN_*.py`, `model_m5.py` (физика поверх ядра) | **A** (ведущий физику) — БЕЗ хэндоффа к B | read; B — ревьюер эквивалентности |
| Геометрия образцов `research/experiments/geometry.py` (вынести `fit_lib.GEOMETRY` сюда) | **A** | read (ORCH, все) |
| `attenuator_app/**`, `docs/attenuator_app/**` | C | read |
| `research/paper/**` | D | read |
| `../THz-WGP-Analysis/litrev/**` (кросс-репо), `research/papers/**` | **L** | read |
| `coordination/**`, `CLAUDE.md`, `research/state.json` (main), `research/RESEARCH_PLAN.md` | ORCH | read; предложения — в ACTIVITY |

## Общие (shared) файлы и правила доступа
| Файл | Правило |
|---|---|
| `research/logs/RESEARCH_LOG.md` | **append-only**, все сессии дописывают с тегом `[<ID>]`; не редактировать чужие записи |
| `research/hypotheses/HYPOTHESES.md` | с 2026-07-28 — ЭКСКЛЮЗИВНАЯ зона A (см. выше) + ORCH; прочим read, правка через HANDOFFS |
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

## Граница B↔A: ПО СЛОЮ, а не по каталогу (решение ORCH, 2026-07-28)
Инвариант, который защищаем, лежит **по слоям**, а не по каталогам:
- **B = общее численное ядро** (`model_core`, `fit_lib`, `verify_*`) + ревьюер эквивалентности.
- **A = скрипты-эксперименты/гипотезы поверх ядра** (`tN_*.py`, `model_m5`), пишет их сам, без хэндоффа.
- **Инвариант CHARTER:** второй копии Бланко/грида/невязки/билдера не заводить — звать `model_core`;
  форк только через хэндофф с обоснованием. Проверяется машинно в `verify_model_core.py` (B добавляет
  поиск прямых вызовов `compute_t_perp`/`compute_theoretical_grid_2d` вне ядра).
- **`GEOMETRY` вынести** из `fit_lib.py` в `research/experiments/geometry.py` (данные образцов, зона A).
- B параметризует `build_experiment(data, angles_limit=None)` (дефолт = `config.ANGLES_LIMIT_2D`),
  A зовёт с `angles_limit=(-90,90)` — вместо шестой копии фильтра углов; повод схлопнуть 5 существующих.

## Конфликтные горячие точки (беречь особо)
- `fit_lib.py`, `model_core.py` — правит только B. A/C/D читают.
- `RESEARCH_LOG.md` — только append, тег сессии.
- `research/state.json` (main) vs `research/two_wgp/state.json` — не путать; каждая сессия свой.
- git `main` — при общем каталоге коммитят по очереди с анонсом; при worktree — мерж через ORCH.
