# BOARD — живой статус сессий

Каждая сессия обновляет СВОЮ строку: при старте, смене задачи, «пульсе» (раз в несколько шагов) и
остановке. Формат времени — ISO (локальное). Статус: `running` / `idle` / `blocked` / `stopped`.

| Сессия | Статус | Текущая задача | План (кратко) | Обновлено (heartbeat) |
|---|---|---|---|---|
| **A** two-WGP | idle | A3_refit_series / A4_calc | фит series1-3 (P16/D11 приор); софт затухания | 2026-07-27 (seed) |
| **B** core/refactor | idle | model_core стабилен | ускорение фитов, дедуп; по запросу | 2026-07-27 (seed) |
| **C** app | idle | attenuator_app прототип | UI/логика калькулятора аттенюатора | 2026-07-27 (seed) |
| **D** paper | idle | DRAFT из SYNTHESIS | Abstract/Methods/Results | 2026-07-27 (seed) |
| **ORCH** control | not_started | — | мониторинг доски, мержи, конфликты | — |

## Легенда владения (полностью — OWNERSHIP.md)
A→`research/two_wgp/` · B→`research/experiments/model_core.py`+shared-code · C→`attenuator_app/` ·
D→`research/paper/` · ORCH→`coordination/`,`CLAUDE.md`,main `research/state.json`.

## Активные «замки» намерений
См. `coordination/sessions/*.md`. Перед правкой общего файла — проверь их.
