# Implementation Plan: Интеграция LMFIT для статистического анализа

> **Для агента IDE**: Следуй этому плану для создания альтернативного конвейера оптимизации на базе `lmfit`. Твоя цель — не просто найти минимум, а извлечь матрицы ковариации, стандартные ошибки и доверительные интервалы для физических параметров.

---

## 1. Концепция и обоснование
Вместо использования стандартного `scipy.optimize.minimize` (Nelder-Mead), который возвращает только точку минимума, мы создадим альтернативный модуль `optimizer_lmfit.py`. 
Пакет `lmfit` использует алгоритм Левенберга-Марквардта (Levenberg-Marquardt), который требует от целевой функции возвращать **вектор невязок**, а не скалярную сумму. 
Склеив взвешенные невязки амплитуды и фазы в один вектор, мы заставим `lmfit` автоматически рассчитать $\chi^2$, стандартные ошибки (standard error) и матрицу корреляций параметров.

## 2. Предлагаемые изменения

### Шаг 1: Создание `unified_optimizer/optimizer_lmfit.py`
Создай новый файл, который будет содержать логику `lmfit`.
```python
import numpy as np
import lmfit
from unified_optimizer import config
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d, find_auto_water_mask

def residual_2d_complex(params, angles_val, analysis_freqs, exp_trans_2d, valid_mask):
    """
    Функция невязки для lmfit. Возвращает 1D массив взвешенных отклонений.
    """
    p_um = params['P_um'].value
    d_um = params['D_um'].value
    loss_factor = params['loss_factor'].value
    gamma = params['gamma'].value
    angle_offset = params['angle_offset'].value
    tau_ps = params['tau_ps'].value

    p = p_um * 1e-6
    d = d_um * 1e-6

    # Штраф, если d >= p (хотя lmfit сам поддерживает bounds, это для надежности)
    if d >= p:
        return np.ones(np.sum(valid_mask)) * 1e6

    theo_complex = compute_theoretical_grid_2d(
        angles_val, analysis_freqs, p, d, loss_factor, angle_offset, tau_ps, gamma=gamma
    )

    exp_masked = exp_trans_2d[valid_mask]
    theo_masked = theo_complex[valid_mask]

    # Амплитудная невязка
    amp_residual = np.abs(exp_masked) - np.abs(theo_masked)
    
    # Фазовая невязка (корректный unwrap разности)
    phase_residual = np.angle(exp_masked) - np.angle(theo_masked)
    phase_residual = np.arctan2(np.sin(phase_residual), np.cos(phase_residual))

    # Веса
    W_AMP = 1.0
    W_PHASE = np.sqrt(0.1) # Корень из веса, так как lmfit сам возведет в квадрат

    # Склеиваем в один 1D массив (LMFIT сам возведет элементы в квадрат и сложит)
    return np.concatenate([amp_residual * W_AMP, phase_residual * W_PHASE])

def run_lmfit_2d(data_dict, angles, dataset_name=""):
    """
    Подготовка данных и запуск LMFIT. Возвращает объект результата (MinimizerResult).
    (Скопировать логику подготовки масок и частот из optimizer_2d.py)
    """
    # ... логика подготовки массивов (как в optimizer_2d.py) ...
    
    # Инициализация параметров LMFIT
    params = lmfit.Parameters()
    p_fixed_um = (config.P_FIXED * 1e6) if config.P_FIXED is not None else 16.0
    
    # Если P фиксирован, устанавливаем vary=False
    params.add('P_um', value=p_fixed_um, vary=(config.P_FIXED is None))
    params.add('D_um', value=5.0, min=1.0, max=p_fixed_um - 0.5)
    params.add('loss_factor', value=0.3, min=0.0, max=5.0)
    params.add('gamma', value=1.0, min=0.1, max=3.0, vary=config.OPTIMIZE_GAMMA)
    params.add('angle_offset', value=0.0, min=-10.0, max=10.0)
    params.add('tau_ps', value=0.0, min=-1.0, max=1.0)
    
    # ... определение valid_mask ...

    mini = lmfit.Minimizer(residual_2d_complex, params, fcn_args=(angles_val, analysis_freqs, exp_trans_2d, valid_mask))
    result = mini.minimize(method='leastsq') # Levenberg-Marquardt
    
    return result
```

### Шаг 2: Создание скрипта `scripts/run_lmfit_analysis.py`
Создай скрипт для запуска LMFIT на одной лучшей серии данных (`356att` или `series3`).
Скрипт должен:
1. Загрузить датасет через `DataManager`.
2. Вызвать `run_lmfit_2d`.
3. Распечатать отчет: `print(lmfit.fit_report(result))`
4. Попробовать рассчитать точные доверительные интервалы: `ci = lmfit.conf_interval(mini, result)` и вывести их с помощью `lmfit.printfuncs.report_ci(ci)`.
5. Сохранить текстовый лог отчета в `results/lmfit_report.txt`.

## 3. Требования к зависимостям
Агенту IDE потребуется убедиться, что `lmfit` установлен в окружении проекта. Если нет, он должен добавить его в `requirements.txt` (если такой есть) или просто выполнить `pip install lmfit`.

## 4. Ожидаемый результат (Верификация)
После прогона `run_lmfit_analysis.py` мы получим файл отчета, в котором будет указано:
- Матрица корреляций (чтобы увидеть, насколько сильно связаны, например, `D_um` и `gamma`).
- Стандартные ошибки (например, $D_{eff} = 4.40 \pm 0.05$ мкм). Это докажет физическую значимость извлеченного эффективного диаметра.
