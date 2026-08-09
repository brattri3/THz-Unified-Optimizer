# Интеграция и воспроизведение оригинальной 2D Nelder-Mead оптимизации

Мы планируем интегрировать поддержку прецизионной 2D спектрально-угловой Nelder-Mead оптимизации со степенным законом затухания и фиксацией шага решетки в единый фреймворк `unified_optimizer`.

## User Review Required

> [!NOTE]
> Наш анализ показал, что свободная оптимизация периода $P$ по ТГц-спектрам всегда упирается в границы из-за вырожденности модели Бланко в субволновом режиме (flat valley phenomenon).
> Паспортные параметры ($P_{\text{eff}} = 15.50$ мкм, $D_{\text{eff}} = 5.67$ мкм) были получены при фиксированном $P = 15.50$ мкм (измеренном независимо по микрофотографиям) и использовании высокостабильной Серии 1.
> Мы добавим в фреймворк возможность фиксации геометрии и степенного закона потерь, чтобы воспроизвести эти паспортные результаты.

## Proposed Changes

### [1. Конфигурация] [config.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/config.py)
Добавим настройки для 2D-оптимизации:
```python
# Настройки прецизионной 2D оптимизации
F_MAX_2D = 1.5              # Ограничение частоты для 2D фитинга (ТГц)
P_FIXED = 15.50e-6          # Фиксированный шаг решетки (м) или None для свободного подбора
USE_POWER_LAW = True        # Использовать степенной закон затухания
OPTIMIZE_GAMMA = True       # Подбирать показатель степени gamma или зафиксировать на GAMMA_DEFAULT
GAMMA_DEFAULT = 1.58        # Показатель степени по умолчанию
ANGLES_LIMIT_2D = (0.0, 90.0) # Диапазон углов для 2D фитинга или None для всех
```

### [2. Оптимизатор] [optimizer_2d.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/optimizer_2d.py)
- Перепишем функцию расчета теоретической сетки `compute_theoretical_grid_2d` с поддержкой степенного закона затухания амплитуды:
  $$t_{\text{eff}} = t \cdot e^{-0.5 \cdot (\alpha_{\text{dB}} / 4.343) \cdot \nu^\gamma}$$
- Перепишем `optimize_2d_spectral` для поддержки:
  - Фильтрации углов по `ANGLES_LIMIT_2D`.
  - Частотного диапазона до `F_MAX_2D`.
  - Фиксации периода `P` на `P_FIXED`.
  - Оптимизации или фиксации `gamma`.
  - Соответствующей размерности вектора подбора в Nelder-Mead (от 2 до 5 параметров).

### [3. CLI-интерфейс] [run_pipeline.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/run_pipeline.py)
- Обновим вывод 2D результатов в консоль, чтобы отображать все найденные параметры (P, D, loss_db, gamma, angle_offset).

### [4. Отчетность] [analytics.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/analytics.py)
- Обновим сравнительную таблицу отчета для корректного отображения параметров степенного закона (`gamma`) и периода.

---

## Verification Plan

### Automated Tests
- Запустим пайплайн для Серии 2 и Серии 1:
  ```powershell
  python -m unified_optimizer.run_pipeline --dataset series1 --method 2d
  python -m unified_optimizer.run_pipeline --dataset series2 --method 2d
  ```
- Проверим сходимость параметров `P_eff`, `D_eff` и `gamma` на Серии 1 к паспортным значениям.
