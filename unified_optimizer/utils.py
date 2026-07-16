import numpy as np
from scipy.signal import medfilt

def normalize_angle(angle: float) -> float:
    return ((angle + 180.0) % 360.0) - 180.0

def find_auto_water_mask(F: np.ndarray, Y: np.ndarray, window_width=0.5, window_step=0.25, threshold_sigma=3.0, cut_width=0.05):
    """
    Автоматический поиск статистических выбросов (линий воды) с использованием бегущего окна.
    Возвращает:
        global_mask: булев массив, где True - хорошие точки, False - выбросы
        masked_intervals: список кортежей (start_f, end_f) вырезанных интервалов
    """
    masked_intervals = []
    global_mask = np.ones(len(F), dtype=bool)
    
    # Детрендинг: вычитание базовой линии с помощью медианного фильтра
    df = np.median(np.diff(F))
    if df > 0:
        kernel_size = int(window_width / df)
        if kernel_size % 2 == 0:
            kernel_size += 1
        # Убедимся, что kernel_size >= 3
        kernel_size = max(3, kernel_size)
        baseline = medfilt(Y, kernel_size)
    else:
        baseline = Y
        
    Y_detrend = Y - baseline
    
    iteration = 0
    while True:
        max_outlier_val = -1.0
        max_outlier_idx = -1
        
        start_f = F[0]
        end_f = F[-1]
        
        curr_f = start_f
        while curr_f < end_f:
            # Маска для текущего окна с учетом уже вырезанных точек
            w_mask = (F >= curr_f) & (F < curr_f + window_width) & global_mask
            if np.sum(w_mask) > 5:
                window_mean = np.mean(Y_detrend[w_mask])
                window_std = np.std(Y_detrend[w_mask])
                
                if window_std > 1e-12:
                    deviations = np.abs(Y_detrend[w_mask] - window_mean)
                    if np.max(deviations) > threshold_sigma * window_std:
                        idx_in_window = np.argmax(deviations)
                        global_idx = np.where(w_mask)[0][idx_in_window]
                        dev_val = deviations[idx_in_window]
                        
                        if dev_val > max_outlier_val:
                            max_outlier_val = dev_val
                            max_outlier_idx = global_idx
            
            curr_f += window_step
            
        if max_outlier_idx != -1:
            outlier_f = F[max_outlier_idx]
            cut_mask = (F >= outlier_f - cut_width) & (F <= outlier_f + cut_width)
            global_mask[cut_mask] = False
            masked_intervals.append((outlier_f - cut_width, outlier_f + cut_width))
            iteration += 1
        else:
            break
            
    # Сортировка и объединение перекрывающихся интервалов
    if not masked_intervals:
        return global_mask, []
        
    masked_intervals.sort(key=lambda x: x[0])
    merged = [masked_intervals[0]]
    for current in masked_intervals[1:]:
        previous = merged[-1]
        if current[0] <= previous[1]:
            merged[-1] = (previous[0], max(previous[1], current[1]))
        else:
            merged.append(current)
            
    return global_mask, merged

def compute_numerical_hessian(func, x0, h=1e-3):
    """
    Вычисляет матрицу Гессиана для функции func в точке x0 численно методом центральных разностей.
    """
    n = len(x0)
    hessian = np.zeros((n, n))
    x0 = np.array(x0, dtype=float)
    
    # Вектор шагов h
    h_vec = np.zeros(n)
    for i in range(n):
        h_vec[i] = h * max(abs(x0[i]), 1.0)
        
    for i in range(n):
        for j in range(n):
            if i == j:
                # Диагональные элементы
                x_plus = x0.copy()
                x_plus[i] += h_vec[i]
                x_minus = x0.copy()
                x_minus[i] -= h_vec[i]
                
                f_plus = func(x_plus)
                f_minus = func(x_minus)
                f_0 = func(x0)
                
                hessian[i, i] = (f_plus - 2 * f_0 + f_minus) / (h_vec[i] ** 2)
            elif i < j:
                # Внедиагональные элементы
                x_pp = x0.copy()
                x_pp[i] += h_vec[i]
                x_pp[j] += h_vec[j]
                
                x_pm = x0.copy()
                x_pm[i] += h_vec[i]
                x_pm[j] -= h_vec[j]
                
                x_mp = x0.copy()
                x_mp[i] -= h_vec[i]
                x_mp[j] += h_vec[j]
                
                x_mm = x0.copy()
                x_mm[i] -= h_vec[i]
                x_mm[j] -= h_vec[j]
                
                f_pp = func(x_pp)
                f_pm = func(x_pm)
                f_mp = func(x_mp)
                f_mm = func(x_mm)
                
                val = (f_pp - f_pm - f_mp + f_mm) / (4 * h_vec[i] * h_vec[j])
                hessian[i, j] = val
                hessian[j, i] = val
                
    return hessian

def estimate_errors_from_hessian(hessian, min_loss, num_points):
    """
    Оценивает стандартные ошибки параметров на основе Гессиана.
    """
    try:
        inv_h = np.linalg.inv(hessian)
        # Ковариационная матрица Sigma = (2 * min_loss / M) * H^-1
        scale = 2.0 * min_loss / max(num_points, 1)
        cov = inv_h * scale
        variances = np.diag(cov).copy()
        variances[variances < 0] = 0.0
        return np.sqrt(variances), cov
    except np.linalg.LinAlgError:
        return np.zeros(len(hessian)), None
