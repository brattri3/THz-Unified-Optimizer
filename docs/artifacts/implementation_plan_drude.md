# Implementation Plan: Drude Model & Diagnostics

> **Для агента IDE**: Строго выполняй шаги по порядку. После каждого шага проводи верификацию. В конце зафиксируй изменения через git.

---

## Шаг 1 — Скрипт диагностики симметрии углов
**Задача**: Создать скрипт, который проверит, насколько данные при углах $\theta$ и $-\theta$ (или $360-\theta$) симметричны для `series1` и `series2`.

**Действие**:
Создай файл `scripts/diagnose_symmetry.py`:
```python
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from unified_optimizer.data_manager import DataManager

def diagnose():
    manager = DataManager()
    for ds in ['series1', 'series2']:
        print(f"\n--- Диагностика {ds} ---")
        try:
            data = manager.get_data_for_dataset(ds)
        except Exception:
            continue
            
        angles = sorted(list(data.keys()))
        positive_angles = [a for a in angles if 0 <= a <= 90]
        
        for pos_a in positive_angles:
            neg_a = -pos_a
            if neg_a in data and pos_a != 0:
                t_s_pos, E_s_pos, _, _ = data[pos_a]
                t_s_neg, E_s_neg, _, _ = data[neg_a]
                
                # Сравниваем интегралы энергии
                int_pos = np.trapezoid(E_s_pos**2, t_s_pos)
                int_neg = np.trapezoid(E_s_neg**2, t_s_neg)
                
                diff_percent = abs(int_pos - int_neg) / ((int_pos + int_neg)/2) * 100
                print(f"Угол {pos_a:2d} vs {neg_a:3d}: разница {diff_percent:.1f}%")

if __name__ == "__main__":
    diagnose()
```
**Верификация**: Выполни `python scripts/diagnose_symmetry.py`. Если разница превышает 10%, данные несимметричны.

---

## Шаг 2 — Внедрение модели Друде в `model_blanco.py`
**Задача**: Заменить идеальный импеданс металла на физичный импеданс вольфрама по модели Друде. Феноменологическое затухание (рассеяние на дефектах) при этом *сохраняется* в оптимизаторе!

**Действие 1**: Открой `unified_optimizer/model_blanco.py`.
Добавь после `C_LIGHT = 3e8` (строка 5) новую константу:
```python
Z0 = 376.7303  # Импеданс вакуума
```

Добавь новую функцию для расчёта нормализованного импеданса Друде:
```python
def get_drude_impedance_normalized(freq_thz: float) -> complex:
    """Нормализованный поверхностный импеданс вольфрама Z_s / Z_0."""
    if freq_thz is None or freq_thz <= 0:
        return 0j
    omega = freq_thz * 2 * np.pi * 1e12
    sigma0 = 1.8e7  # См/м (Вольфрам)
    tau = 8.0e-15   # с
    sigma_omega = sigma0 / (1.0 - 1j * omega * tau)
    mu0 = 4 * np.pi * 1e-7
    Z_s = np.sqrt(1j * omega * mu0 / sigma_omega)
    return Z_s / Z0
```

**Действие 2**: Обнови `compute_t_perp` и `compute_t_par` для приема `freq_thz` и добавления импеданса.
Измени сигнатуры и код:
```python
def compute_t_perp(p_over_lambda: float, d_over_p: float, N: int = 15, freq_thz: float = None) -> complex:
    """Амплитудный коэффициент пропускания перпендикулярной поляризации."""
    fa = compute_fa(p_over_lambda, d_over_p, N)
    fb = compute_fb(p_over_lambda, d_over_p, N)
    Z_drude = get_drude_impedance_normalized(freq_thz)
    
    Za = -1j / fa + Z_drude if abs(fa) > EPS else -1e9j
    Zb = -1j / fb + Z_drude if abs(fb) > EPS else -1e9j
    # ... остальной код без изменений ...

def compute_t_par(p_over_lambda: float, d_over_p: float, N: int = 15, freq_thz: float = None) -> complex:
    """Амплитудный коэффициент пропускания параллельной поляризации."""
    fc = compute_fc(p_over_lambda, d_over_p, N)
    fd = compute_fd(p_over_lambda, d_over_p)
    Z_drude = get_drude_impedance_normalized(freq_thz)
    
    Zc = 1j * fc + Z_drude
    Zd = -1j * fd + Z_drude
    # ... остальной код без изменений ...
```

**Действие 3**: В `unified_optimizer/optimizer_2d.py` (строка ~40), передай частоту `f` внутрь `compute_t_perp` и `compute_t_par`:
**Было**:
```python
        t_perp_arr.append(model_blanco.compute_t_perp(p_over_lambda, d_over_p, N))
        t_par_arr.append(model_blanco.compute_t_par(p_over_lambda, d_over_p, N))
```
**Стало**:
```python
        t_perp_arr.append(model_blanco.compute_t_perp(p_over_lambda, d_over_p, N, freq_thz=f))
        t_par_arr.append(model_blanco.compute_t_par(p_over_lambda, d_over_p, N, freq_thz=f))
```
**Внимание**: Феноменологическая функция потерь (с `loss_factor` и `gamma`) в `optimizer_2d.py` остаётся *как есть*! Она теперь моделирует только диффузное рассеяние на дефектах.

---

## Шаг 3 — Корректировка Global Average
**Задача**: Исключить сомнительные серии из глобального усреднения.
**Файл**: `scripts/run_overnight_pipeline.py`

Найди блок `if __name__ == "__main__":` в самом низу.
**Было**:
```python
    manager = DataManager()
    all_datasets = manager.get_available_datasets()
```
**Стало**:
```python
    manager = DataManager()
    all_datasets = manager.get_available_datasets()
    
    # Оставляем только физически надежные серии для Global Average
    reliable_datasets = ['356att', 'series3']
    all_datasets = [ds for ds in all_datasets if ds in reliable_datasets]
```

**Верификация**: Запусти `python scripts/run_overnight_pipeline.py --dry-run`. Убедись, что выполняются только `356att` и `series3`.

---

## Финал
1. Запусти `python scripts/run_overnight_pipeline.py` без флага `--dry-run`.
2. Добавь изменения в `CHANGELOG.md` (Внедрена модель Друде, добавлен скрипт диагностики).
3. Сделай `git add .`, `git commit -m "Внедрена модель Друде для вольфрама"`, `git push`.
