# Implementation Plan: THz-Unified-Optimizer Fixes

> **Для агента**: Прочти `CHANGELOG.md`, выполни `git pull`, затем строго выполняй задачи по порядку. После каждой задачи — верификация. В конце — `git commit` и `git push`.

---

## Задача 1 — Исключить series4 и series5 из пайплайна
**Файл**: `scripts/run_overnight_pipeline.py`, строки 60–63

**Было** (строки 60–63):
```python
        angles = sorted(list(data_dict.keys()))
        logging.info(f"Загружено точек по углам: {len(angles)}")

        results[ds] = {}
```

**Стало**:
```python
        angles = sorted(list(data_dict.keys()))
        logging.info(f"Загружено точек по углам: {len(angles)}")

        MIN_ANGLES = 5  # Минимально необходимое число угловых точек
        if len(angles) < MIN_ANGLES:
            logging.warning(f"[{ds}] ПРОПУСК: только {len(angles)} угловых точек (требуется >= {MIN_ANGLES}). Серия исключена.")
            continue

        results[ds] = {}
```

**Верификация**:
```bash
python scripts/run_overnight_pipeline.py --dry-run
```
В логе должны появиться строки `ПРОПУСК` для series4 и series5.

---

## Задача 2 — Исправить верхнюю границу параметра D в 2D-оптимизаторе
**Файл**: `unified_optimizer/optimizer_2d.py`, строка 118

**Было** (строка 118):
```python
    free_params.append(('D_um', config.D_DEFAULT * 1e6, (1.0, 25.0)))
```

**Стало**:
```python
    p_fixed_um = (config.P_FIXED * 1e6) if config.P_FIXED is not None else 40.0
    free_params.append(('D_um', config.D_DEFAULT * 1e6, (1.0, p_fixed_um - 0.5)))
```

Добавить эту строку **перед** строкой 118 (или заменить строку 118 на две):
- `p_fixed_um` вычисляется на основе `config.P_FIXED` (сейчас = 15.50 мкм)
- Теперь D ограничен сверху значением 15.00 мкм — вырожденное решение D=P невозможно

**Верификация**: после правки убедиться, что строка 148 (`if d_um >= p_um`) больше не является единственным барьером:
```python
# Проверить что bounds сами не пропускают D >= P_FIXED
# Строка 148 остаётся как дополнительная защита — не удалять
```

---

## Задача 3 — Добавить логирование сходимости (ConvergenceLogger)

### 3a. В файл `unified_optimizer/optimizer_2d.py`

**Вставить** класс и его использование. В начало функции `optimize_2d_spectral` (после строки 71, перед строкой 72):

В **начало файла** (после `import time`, строка 3) добавить:
```python
import logging
```

Перед строкой `res = minimize(` (строка 180) вставить:

```python
    # --- Логирование сходимости ---
    _convergence_log = []
    _iter_counter = [0]

    def _callback_2d(xk):
        _iter_counter[0] += 1
        if _iter_counter[0] % 100 == 0:
            loss_val = loss_function_2d(xk)
            params_str = ', '.join(
                f"{name}={val:.4f}" for val, (name, _, _) in zip(xk, free_params)
            )
            logging.info(f"  [2D iter={_iter_counter[0]:4d}] loss={loss_val:.6f} | {params_str}")
        _convergence_log.append(_iter_counter[0])

```

Заменить строки 180–185 (`res = minimize(...)`) на:
```python
    res = minimize(
        loss_function_2d,
        initial_guess,
        method='Nelder-Mead',
        callback=_callback_2d,
        options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5}
    )
    logging.info(f"  [2D] Завершено: success={res.success}, nit={res.nit}, nfev={res.nfev}")
    logging.info(f"  [2D] Причина: {res.message}")
```

Добавить проверку границ **перед `return`** (перед строкой `return {`, строка 203):
```python
    # Предупреждение при попадании на границы
    for val, (name, _, bounds) in zip(res.x, free_params):
        if abs(val - bounds[0]) < 0.01 or abs(val - bounds[1]) < 0.01:
            logging.warning(f"  [2D] ВНИМАНИЕ: параметр '{name}'={val:.4f} на границе {bounds}!")
```

### 3b. В файл `unified_optimizer/optimizer_1d.py`

В **начало файла** добавить:
```python
import logging
```

После строки 70 (`res1 = minimize(...)`) добавить:
```python
    logging.info(f"  [1D stage1] fun={res1.fun:.4f}, nit={res1.nit}, x={res1.x.round(3).tolist()}")
```

После строки 82 (`res2 = minimize(...)`) добавить:
```python
    logging.info(f"  [1D stage2] fun={res2.fun:.4f}, nit={res2.nit}, x={res2.x.round(3).tolist()}")
```

После строки 94 (`res3 = minimize(...)`) добавить:
```python
    logging.info(f"  [1D stage3] fun={res3.fun:.4f}, nit={res3.nit}, x={res3.x.round(3).tolist()}")
    # Предупреждение при попадании на границы
    names3 = ['P_um', 'D_um', 'alpha', 'gamma', 'theta_offset', 'eps_floor']
    for name, val, bounds in zip(names3, res3.x, b3):
        if abs(val - bounds[0]) < 0.01 or abs(val - bounds[1]) < 0.01:
            logging.warning(f"  [1D] ВНИМАНИЕ: параметр '{name}'={val:.4f} на границе {bounds}!")
```

**Верификация**:
```bash
python scripts/run_overnight_pipeline.py --dry-run
```
В логе должны появиться строки `[2D iter=...]` и `[1D stage1] fun=...`.

---

## Задача 4 — Разделить амплитудную и фазовую невязку в 2D целевой функции
**Файл**: `unified_optimizer/optimizer_2d.py`, строки 165–178

**Было** (строки 165–178):
```python
        diff_complex = exp_trans_2d[valid_mask] - theo_complex[valid_mask]
        diff_abs = np.abs(diff_complex)
        
        # Удаление статистических выбросов (> 3 сигма) для устойчивости фиттинга
        mean_diff = np.mean(diff_abs)
        std_diff = np.std(diff_abs)
        inliers = diff_abs < (mean_diff + 3 * std_diff)
        
        if np.sum(inliers) == 0:
            return 1e6
            
        rmse_complex = np.sqrt(np.mean(diff_abs[inliers]**2))
        
        return rmse_complex
```

**Стало**:
```python
        exp_masked = exp_trans_2d[valid_mask]
        theo_masked = theo_complex[valid_mask]

        # Амплитудная невязка
        amp_residual = np.abs(exp_masked) - np.abs(theo_masked)

        # Фазовая невязка (корректный unwrap разности)
        phase_residual = np.angle(exp_masked) - np.angle(theo_masked)
        phase_residual = np.arctan2(np.sin(phase_residual), np.cos(phase_residual))

        # 3-сигма отсечение по амплитудной невязке
        amp_abs = np.abs(amp_residual)
        mean_a, std_a = np.mean(amp_abs), np.std(amp_abs)
        inliers = amp_abs < (mean_a + 3 * std_a)

        if np.sum(inliers) == 0:
            return 1e6

        # Раздельная функция потерь: амплитуда + фаза с весами
        W_AMP = 1.0    # вес амплитудной невязки
        W_PHASE = 0.1  # вес фазовой невязки (меньше — фаза шумнее)
        loss_amp = np.mean(amp_residual[inliers]**2)
        loss_phase = np.mean(phase_residual[inliers]**2)
        return W_AMP * loss_amp + W_PHASE * loss_phase
```

**Верификация**:
```bash
python scripts/run_overnight_pipeline.py --dry-run
```
Прогон должен завершиться без ошибок. Значение `fun` в логе теперь будет другого масштаба (не RMSE, а взвешенная сумма MSE) — это нормально.

---

## Финальная верификация

Запустить полный прогон на двух лучших сериях:
```bash
python scripts/run_overnight_pipeline.py
```

Проверить в `results/overnight_results.json`:
- **series4 и series5 отсутствуют** в JSON
- **D для 356att и series3** не равно 15.50 мкм (не вырождено)
- **D < 15.00 мкм** для всех серий (граница сработала)

Проверить в `results/overnight_execution.log`:
- Есть строки `[2D iter=...]` — логирование работает
- Есть строки `[1D stage1] fun=...` — логирование работает
- Нет строк `ВНИМАНИЕ` для D параметра — вырождение устранено

---

## После завершения всех задач

1. Обновить `CHANGELOG.md` — добавить запись о внесённых правках
2. Скопировать этот файл в `docs/artifacts/` (он уже там: `implementation_plan_fixes.md`)
3. Выполнить:
```bash
git add .
git commit -m "Исправление 4 критических багов: bounds D, ConvergenceLogger, разделение amp/phase невязки, фильтр series4/5"
git push
```
