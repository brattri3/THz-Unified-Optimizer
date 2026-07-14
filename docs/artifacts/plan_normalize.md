# План интеграции нормализации углов решетки поляризатора

Этот план описывает создание вспомогательного модуля `utils.py` с функцией нормализации углов и интеграцию этой функции в менеджер данных, чтобы привести углы к диапазону от -180° до 180°.

## User Review Required

> [!IMPORTANT]
> Мы создадим новый файл `unified_optimizer/utils.py` с функцией `normalize_angle`. Эта функция будет интегрирована в функцию `parse_filename` в `data_manager.py`. Все углы в пайплайне (например, 270° -> -90°) будут автоматически приведены к диапазону [-180°, 180°].

## Proposed Changes

### 1. Создание модуля утилит

#### [NEW] [utils.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/utils.py)
Создание файла со следующим кодом:
```python
import numpy as np

def normalize_angle(angle: float) -> float:
    return ((angle + 180.0) % 360.0) - 180.0
```

### 2. Интеграция нормализации углов

#### [MODIFY] [data_manager.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/data_manager.py)
- Импорт `normalize_angle` из `unified_optimizer.utils`.
- Применение `normalize_angle` к распарсенному углу `angle_deg` в функции `parse_filename`.

### 3. Обновление автотестов

#### [MODIFY] [test_data_manager.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/tests/test_data_manager.py)
- Добавление теста на проверку нормализации углов, выходящих за границы [-180, 180] (например, 270 -> -90).

## Verification Plan

### Automated Tests
- Запуск тестов:
  ```powershell
  python -m unittest discover -s unified_optimizer/tests
  ```

### Manual Verification
- Запуск пайплайна для датасета `series1` (в котором углы 270°-350° должны превратиться в -90°...-10°):
  ```powershell
  python -m unified_optimizer.run_pipeline --dataset series1 --method 1d
  ```
- Проверка вывода загруженных углов на консоль (углы больше 180 градусов должны быть заменены на отрицательные).
