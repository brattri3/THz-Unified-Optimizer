# CHANGELOG — рефактор `research/experiments` в кэширующее ядро `model_core`

Ветка: `main`. Гардрейлы `research/RESEARCH_PLAN.md §6` соблюдены:
`unified_optimizer/` и `data_pool/` — **только чтение**, все записи под `research/`,
никаких `rm -rf` / `git push` / сетевых операций.

---

## 2026-07-27 — Итерация 1: ядро `model_core` + переключение t6 / t8 / model_m5 / t20

### Что было
Каждый фит-скрипт нёс свою почти дословную копию `compute_theoretical_grid_2d`,
и каждая копия на **каждом `nfev`** заново крутила питоновский цикл по частотам
с вызовами `model_blanco.compute_t_perp` / `compute_t_par` (внутри — ещё цикл
`m=1..15`). При этом Бланко зависит только от `(сетка ν, P, D, N, drude)` и
**не** зависит от `loss_factor / gamma / angle_offset / tau_ps / tau_par_ps /
eta0 / eta_exp / delta0 / tau_leak` — то есть при численном якобиане lmfit
пересчитывался идентично для каждого столбца, не трогающего `D`.

Дубли грида: `optimizer_2d.compute_theoretical_grid_2d`,
`t6_hn2.compute_grid_hn2`, `t6_hn8.compute_grid_hn8` — три копии одной формулы.
Дубли невязки: `fit_lib.residual`, `t6_hn2.residual_hn2`, `t6_hn8.residual`,
`model_m5.residual`, `t20_weighted_fit.make_resid` — пять копий.

### Добавлено
**`research/experiments/model_core.py`** — общее ядро:
- `blanco_t(freqs, p, d, N, use_drude)` — **векторизованный по частоте** Бланко
  с LRU-мемоизацией (ключ — точные биты `float`, без округления: lmfit
  возмущает `D` на ~1e-8 относительных, округление ключа склеило бы соседние
  точки и испортило производную). Ёмкость 512, вытеснение LRU.
- `compute_grid(...)` — **единственная** функция грида, строгий супермножество:
  `eta0=0` → `compute_theoretical_grid_2d`; `tau_leak=0` → `compute_grid_hn2`;
  общий случай → `compute_grid_hn8`.
- `grid_from_params(params, ...)` — грид прямо из `lmfit.Parameters`
  (отсутствующие имена → нейтральные значения по умолчанию).
- `complex_residual(...)` / `residual_from_params(...)` — единая
  amplitude+phase невязка, с опциональными по-точечными весами `w`
  (взвешенный вариант t20).
- `cache_stats()` / `clear_cache()`.

**`research/experiments/verify_model_core.py`** — регрессионный сторож.
Держит **замороженную копию** до-рефакторного `compute_grid_hn8`, чтобы
эквивалентность оставалась проверяемой после того, как `t6_hn2`/`t6_hn8` стали
тонкими алиасами. Допуск — **0** (точное равенство `float`).

**`research/experiments/refactor_ab_check.py`** — A/B-харнесс: гоняет 29 фитов
через публичные точки входа (`t6_hn2.fit_variant`, `t6_hn8.fit`,
`model_m5.fit_variant`, взвешенный фит t20) и сравнивает AIC/BIC/redchi/nfev и
все подогнанные параметры до и после. Артефакты — `research/results/refactor/`.

### Изменено (имена `tN_` НЕ переименованы, публичные функции сохранены)
| Файл | Что стало |
|---|---|
| `t6_hn2.py` | `compute_grid_hn2` → алиас `model_core.compute_grid(tau_leak=0)`; `residual_hn2` → `model_core.residual_from_params`; пост-фит грид → `grid_from_params`. Имя `compute_grid_hn2` сохранено — его импортирует `t6_hn6`. |
| `t6_hn8.py` | `compute_grid_hn8` → алиас `model_core.compute_grid`; `residual` → `model_core.residual_from_params`. Имя сохранено — импортируют `model_m5` и `t20`. |
| `model_m5.py` | `residual` → `model_core.residual_from_params`; пост-фит грид → `grid_from_params`. |
| `t20_weighted_fit.py` | `make_resid` → `model_core.residual_from_params(w=...)`. |
| `t8_consolidate.py` | Без правок — вызывает `model_m5`, ускорение получил транзитивно. |
| `t6_hn6.py` | Без правок — вызывает `t6_hn2.compute_grid_hn2`, ускорение транзитивно. |

### Численный контракт: побитово, а не «в пределах допуска»
Векторизация выполняет **те же операции float64 в том же порядке**: суммы по
`m` по-прежнему накапливаются питоновским циклом (порядок суммирования
сохранён), каждый скалярный `if`-guard стал `np.where` по тому же предикату.
Две ловушки, найденные при доведении до точного равенства:

1. **`_cpython_cdiv`** — `get_drude_impedance_normalized` в оригинале считает на
   *Python-complex*, а CPython делит комплексные по Смиту с финальным делением,
   тогда как NumPy умножает на обратное к знаменателю → расхождение ~1 ULP.
   Реализовано векторное деление, побитово совпадающее с CPython.
2. **`_cmul` / `_csq`** — ufunc комплексного умножения NumPy сворачивает
   `ar*br - ai*bi` в FMA на x86-64, а скалярный путь `np.complex128` — нет:
   расхождение на ~45% случайных пар (8730/20000). Разложение на четыре
   отдельных вещественных ufunc-вызова воспроизводит скалярный результат точно
   (0/20000 расхождений). Комплексное **деление**, сложение и `real*complex`
   у скаляра и массива совпадают — им спец-обработка не нужна.

`verify_model_core.py`: **ALL CHECKS EXACT** (t_perp/t_par на 3 диаметрах ×
2 режимах Друде × 2 набора данных; грид против всех трёх легаси-форм включая
`eta0=0`; невязка обычная и взвешенная).

Дополнительно убраны импорты `from t6_hn8 import compute_grid_hn8` в `model_m5`
и `t20` — они обнулились именно этим рефактором. Уже мёртвый до него
`from scipy import stats` (в `t6_hn2`, `model_m5`) намеренно оставлен: вне
области задачи, чтобы не разводить шум в диффе.

### Результат A/B (`refactor_ab_check.py --diff pre post`, повтор `--diff pre post2` после чистки импортов)
```
identical (bit-for-bit): 482
differing              : 0
missing on one side    : 0
AIC entries: 32, differing: 0
```
Все 29 фитов, все 482 числовых поля (AIC, BIC, redchi, chisqr, nfev, nvarys,
ndata и каждый подогнанный параметр) — **побитово идентичны**.

### Ускорение
| Замер | Было | Стало |
|---|---|---|
| Один вызов грида, `D` меняется | 27–45 мс | 0.5–0.7 мс (**46–71×**) |
| Один вызов грида, `D` в кэше | 27–45 мс | 0.04 мс (**620–1015×**) |
| A/B-набор из 29 фитов, end-to-end | 174.7 с | 11.4 с (**15.3×**) |
| `t6_hn8/356att/hn8/fixD=False` | 19.30 с | 0.47 с (**41×**) |
| `t20/356att` взвешенный | 33.28 с | 0.90 с (**37×**) |

Попадания кэша внутри одного фита (356att):
- `M3ref` (`D` свободен): 65 вызовов грида → 21 пересчёт Бланко, **67.7% хитов**;
- `M5_physD_leak` (`D` приколочен): 162 вызова → **1** пересчёт Бланко, 99.4% хитов.

### Побочная находка (НЕ следствие рефактора)
Закоммиченный `research/results/HN8/hn8_ablation.json` **устарел**: он содержит
`P_um: 40.0` для `test_grid_40_20`, тогда как `fit_lib.GEOMETRY` уже переведён
на микрографную GT `P=38.8` (см. комментарий в `fit_lib.py`: «Legacy T1–T14
used nominal 40/20 → re-run those to compare if needed»). Перегенерация даёт
`P_um: 38.8` и AIC `-22945.12 → -22678.91` для `m3/fixD=True`.
Что это не рефактор — доказывает `pre.json`, снятый **до** правок: там уже
`-22678.906527`. Артефакт откачен к `HEAD`; перегенерацию оставляю владельцу.
`research/results/consolidated/m5_variants.json` после перезапуска `t8` —
байт-в-байт прежний (git не видит изменений).

### Не тронуто (осознанно)
- `unified_optimizer/optimizer_2d.py`, `model_blanco.py` — оригиналы, §6.
- `fit_lib.py` — всё ещё зовёт `compute_theoretical_grid_2d` напрямую; на нём
  сидят `t1`, `t9`, `t10`–`t19`. Перевод на `model_core` даст им тот же 40×+,
  но это отдельная итерация со своим A/B (в задании названы t6/t8/model_m5/t20).
- `research/two_wgp/model_2wgp.py` — там своя FIFO-обёртка `blanco_t` со
  скалярным циклом; параллельная сессия правит этот каталог, не трогаю во
  избежание конфликта. Кандидат на перевод на `model_core` позже (даст ~50×).

### Как проверить
```powershell
.venv\Scripts\python.exe research\experiments\verify_model_core.py          # ALL CHECKS EXACT
.venv\Scripts\python.exe research\experiments\refactor_ab_check.py --diff pre post2
```
`pre.json` — снимок до рефактора, `post.json` — после переключения скриптов,
`post2.json` — после чистки мёртвых импортов (финальное состояние). Оба
`post` дают 0 расхождений с `pre`.

---

## 2026-07-27 — Итерация 2: `fit_lib` на ядре

### Ревизия: кто на чём сидит
Аудит показал, что после итерации 1 семейство `t9`–`t19` ускорилось **само**:
`t9`, `t10`, `t11`, `t12`, `t13`, `t16`, `t18` импортируют
`t6_hn8.compute_grid_hn8`, а `t15b` — `model_m5.fit_variant`; и то и другое уже
делегирует в `model_core`. `t17` и `t6_hn7` вообще не строят грид Бланко.

Единственным оставшимся прямым потребителем `compute_theoretical_grid_2d`
оказался `fit_lib` — точка входа `t1_baseline` и `t2_residual_maps`.

### Изменено
| Файл | Что стало |
|---|---|
| `fit_lib.py` | `residual` → `model_core.compute_grid` + `model_core.complex_residual`; пост-фит грид в `fit_model` → `model_core.compute_grid`; импорт `compute_theoretical_grid_2d` убран (остался `get_transmission_spectra`). |

Историческая деталь сохранена намеренно: в `fit_lib.residual` запасное значение
для отсутствующего `gamma` — **2.0**, а не `1.0` из `model_core._GRID_DEFAULTS`,
поэтому `gamma` и `tau_par_ps` передаются явно, а не через `grid_from_params`.
Обе ветки (параметр есть / параметра нет) покрыты тестом.

### Проверка
`verify_model_core.py` дополнен четырьмя случаями
`fit_lib.residual` (drude вкл/выкл × `gamma` как параметр / через запасное 2.0)
против замороженной легаси-формы — все **EXACT**.

A/B (`refactor_ab_check.py --set fitlib`, лестница аблаций M0…M4 на `356att`,
`test_grid_40_20`, `specac` — 15 фитов). Чтобы молчаливая правка пост-фит грида
не проскочила, в снимок добавлены контрольные суммы карты остатков
(`_ampres_sum`, `_ampres_sumsq`, `_theo_absum`):
```
=== diff pre_fitlib -> post_fitlib ===
identical (bit-for-bit): 543
differing              : 0
AIC entries: 15, differing: 0
```
Контрольный повтор набора итерации 1 (`--diff pre post3`): снова 482/482
побитово, 32 AIC без расхождений — правки `fit_lib` не задели ранее
переключённые скрипты. `python -m compileall research/experiments` — OK.

### Ускорение
| Замер | Было | Стало |
|---|---|---|
| Набор M0…M4 × 3 набора данных (15 фитов) | 25.4 с | 3.5 с (**7.3×**) |
| `fit_lib/356att/M4` | 2.15 с | 0.23 с (**9.3×**) |

Попадания кэша внутри фита `fit_lib` (356att): `M0` — 29 вызовов грида на
13 пересчётов Бланко (55.2% хитов), `M4` — 68 на 19 (72.1%). Доля хитов ниже,
чем у `M5_physD_leak` из итерации 1 (99.4%), потому что здесь `D` свободен и
меняется на каждом шаге оптимизатора; выигрыш идёт в основном от векторизации.

### Итог по задаче
Названная область (`model_core` + `t6`/`t8`/`model_m5`/`t20`) закрыта в
итерации 1; итерация 2 добила последний дубль грида в `research/experiments`.
Прямых вызовов `compute_theoretical_grid_2d` в `research/experiments/` больше
нет — остался ровно один грид (`model_core.compute_grid`) и одна невязка
(`model_core.complex_residual`).

### Что осталось за скобками
- `research/two_wgp/model_2wgp.py` — своя FIFO-обёртка со скалярным циклом;
  каталог правит параллельная сессия, не трогаю. Перевод на `model_core` дал бы
  те же ~50×.
- `unified_optimizer/` — оригиналы, §6.
