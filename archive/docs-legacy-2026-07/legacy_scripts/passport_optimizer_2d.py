import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pathlib import Path
import time
import sys

# Добавляем путь к src, если запускаем напрямую
sys.path.append(str(Path(__file__).resolve().parent))

import config
import data_loader
import spectrum
import theoretical

def load_2d_experimental_data(repetition=2, freq_start=0.2, freq_end=1.5):
    """
    Загружает экспериментальные данные для конкретной серии (repetition)
    и рассчитывает спектральный шумовой пол по высокочастотному хвосту фона.
    """
    # 1. Загрузка датасета во временное хранилище
    data_store = data_loader.load_dataset_to_store(config.DATA_DIR)
    
    # 2. Выделяем только нужную серию (по умолчанию Серия 2 как самая чистая)
    keys = sorted([k for k in data_store.keys() if k[0] == 'signal_raw' and k[3] == repetition])
    if not keys:
        raise ValueError(f"В базе данных не найдены сигналы для Серии {repetition}!")
        
    individual_results = {}
    bg_spectra = []
    freqs_common = None
    
    for key in keys:
        _, angle1, angle2, rep = key
        sig_key = ('signal_raw', angle1, angle2, rep)
        bg_key = ('bg_raw', angle1, angle2, rep)
        
        if bg_key in data_store:
            t_sig, E_sig = data_store[sig_key]
            t_bg, E_bg = data_store[bg_key]
            
            freqs, spec_sig, spec_bg, transmission = spectrum.calculate_transmission(
                t_sig, E_sig, t_bg, E_bg
            )
            individual_results[(angle1, angle2)] = transmission
            bg_spectra.append(spec_bg)
            if freqs_common is None:
                freqs_common = freqs
                
    # 3. Расчет спектрального Noise Floor по высокочастотному хвосту фона (выше 3.0 ТГц)
    bg_avg = np.mean(bg_spectra, axis=0)
    tail_indices = np.where(freqs_common >= 3.0)[0]
    if len(tail_indices) == 0:
        tail_indices = np.where(freqs_common >= 2.5)[0]
        
    noise_floor_amplitude = np.mean(bg_avg[tail_indices])
    
    # Спектральный шумовой пол по мощности
    dr_power = (bg_avg / noise_floor_amplitude)**2
    t_noise_floor = 1.0 / dr_power
    
    # 4. Фильтрация частот по рабочему диапазону
    target_indices = np.where((freqs_common >= freq_start) & (freqs_common <= freq_end))[0]
    analysis_freqs = freqs_common[target_indices]
    fit_t_noise = t_noise_floor[target_indices]
    
    # 5. Сборка 2D-матрицы
    angles = sorted(list(individual_results.keys()))
    angles_val = np.array([a[0] for a in angles])
    
    exp_trans_2d = np.zeros((len(angles_val), len(analysis_freqs)))
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]
        
    # Защита от отрицательных и нулевых значений
    exp_trans_db_2d = 10 * np.log10(np.maximum(exp_trans_2d, 1e-12))
    
    return angles_val, analysis_freqs, exp_trans_2d, exp_trans_db_2d, fit_t_noise

def compute_theoretical_grid_2d(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, t_noise, N=15):
    """
    Быстрый аналитический векторизованный расчет 2D-сетки пропускания по Бланко
    с учетом сдвига угла и шумового пола прибора.
    """
    c_light = 3e8
    d_over_p = d / p
    
    # Расчет комплексных коэффициентов t_perp и t_par для каждой частоты
    t_perp_arr = []
    t_par_arr = []
    for f in freqs_thz:
        lambda_m = c_light / (f * 1e12)
        p_over_lambda = p / lambda_m
        t_perp_arr.append(theoretical.compute_t_perp(p_over_lambda, d_over_p, N))
        t_par_arr.append(theoretical.compute_t_par(p_over_lambda, d_over_p, N))
        
    t_perp_arr = np.array(t_perp_arr)
    t_par_arr = np.array(t_par_arr)
    
    # Преобразование углов в радианы с учетом углового смещения нуля шкалы
    adjusted_angles_rad = np.deg2rad(angles_deg - angle_offset)
    
    # Подготовка размерностей для векторизованного произведения
    cos_a = np.cos(adjusted_angles_rad)[:, np.newaxis]
    sin_a = np.sin(adjusted_angles_rad)[:, np.newaxis]
    
    t_perp = t_perp_arr[np.newaxis, :]
    t_par = t_par_arr[np.newaxis, :]
    
    # Применение омических потерь (напряженность поля)
    t_perp_eff = t_perp * np.exp(-0.5 * loss_factor * freqs_thz)[np.newaxis, :]
    t_par_eff = t_par * np.exp(-0.5 * loss_factor * freqs_thz)[np.newaxis, :]
    
    # Амплитуда горизонтальной компоненты на приемнике
    E_out = cos_a**2 * t_perp_eff + sin_a**2 * t_par_eff
    T_ideal = np.abs(E_out)**2
    
    # Добавляем спектральный шумовой пол (по мощности)
    T_mod = T_ideal + t_noise[np.newaxis, :]
    
    return np.clip(T_mod, 0.0, 20.0)

def loss_function_2d(params_scaled, angles, freqs, exp_linear, exp_db, t_noise, w_linear=10.0, w_db=1.0):
    """
    Прецизионная целевая функция для 2D-аппроксимации с масштабированными параметрами.
    Подбирает: [P_um, D_um, loss_factor, angle_offset]
    P_um и D_um передаются в микрометрах (масштаб ~10) для улучшения сходимости Nelder-Mead.
    """
    p_um, d_um, loss_factor, angle_offset = params_scaled
    
    # Физические ограничения параметров
    if p_um <= 5.0 or p_um >= 40.0:
        return 1e8
    if d_um <= 1.0 or d_um >= 25.0:
        return 1e8
    if d_um >= p_um:
        return 1e8 + (d_um - p_um) * 1e9
    if loss_factor < 0.0 or loss_factor > 5.0:
        return 1e8
    if angle_offset < -10.0 or angle_offset > 10.0:
        return 1e8
        
    p = p_um * 1e-6
    d = d_um * 1e-6
        
    theo_linear = compute_theoretical_grid_2d(angles, freqs, p, d, loss_factor, angle_offset, t_noise)
    theo_db = 10 * np.log10(np.maximum(theo_linear, 1e-12))
    
    # Маска достоверных точек (сигнал выше порога шума хотя бы в 1.5 раза)
    valid_mask = exp_linear > 1.5 * t_noise[np.newaxis, :]
    if np.sum(valid_mask) == 0:
        return 1e8
        
    # Расчет ошибок только по достоверным точкам
    diff_lin = exp_linear[valid_mask] - theo_linear[valid_mask]
    rmse_lin = np.sqrt(np.mean(diff_lin**2))
    
    rmse_db = np.sqrt(np.mean((exp_db[valid_mask] - theo_db[valid_mask])**2))
    
    return w_linear * rmse_lin + w_db * rmse_db

def run_2d_optimization(angles, freqs, exp_linear, exp_db, t_noise, w_linear=10.0, w_db=1.0):
    """
    Выполняет Nelder-Mead оптимизацию 4 параметров (P, D, loss_factor, angle_offset)
    в удобных для сходимости масштабах.
    """
    # Начальное приближение в микрометрах: [P_um, D_um, loss_factor, angle_offset]
    initial_guess = [config.P_DEFAULT * 1e6, config.D_DEFAULT * 1e6, config.LOSS_FACTOR_DEFAULT, 0.0]
    bounds = [
        (15.5, 16.5),        # Шаг P (мкм) зажат вокруг паспортных 16.0 мкм
        (1.0, 25.0),        # Диаметр D (мкм)
        (0.0, 5.0),         # Потери loss_factor
        (-5.0, 5.0)         # Сдвиг нуля шкалы углов (град)
    ]
    
    t0 = time.time()
    res = minimize(
        loss_function_2d, 
        initial_guess, 
        args=(angles, freqs, exp_linear, exp_db, t_noise, w_linear, w_db),
        method='Nelder-Mead',
        bounds=bounds,
        options={'maxiter': 500, 'xatol': 1e-5, 'fatol': 1e-5}
    )
    t1 = time.time()
    
    opt_p_um, opt_d_um, opt_loss, opt_offset = res.x
    opt_p = opt_p_um * 1e-6
    opt_d = opt_d_um * 1e-6
    
    # Финальные погрешности
    theo_linear = compute_theoretical_grid_2d(angles, freqs, opt_p, opt_d, opt_loss, opt_offset, t_noise)
    theo_db = 10 * np.log10(np.maximum(theo_linear, 1e-12))
    
    valid_mask = exp_linear > 1.5 * t_noise[np.newaxis, :]
    rmse_lin = np.sqrt(np.mean((exp_linear[valid_mask] - theo_linear[valid_mask])**2))
    rmse_db = np.sqrt(np.mean((exp_db[valid_mask] - theo_db[valid_mask])**2))
    
    info = {
        'success': res.success,
        'message': res.message,
        'fun': res.fun,
        'nfev': res.nfev,
        'time_sec': t1 - t0,
        'p_um': opt_p_um,
        'd_um': opt_d_um,
        'loss_factor': opt_loss,
        'angle_offset': opt_offset,
        'rmse_lin_pct': rmse_lin * 100,
        'rmse_db': rmse_db
    }
    return [opt_p, opt_d, opt_loss, opt_offset], info

def plot_2d_results(angles, freqs, exp_linear, exp_db, t_noise, opt_params, title_suffix=""):
    """
    Строит и сохраняет 2D-карты результатов фиттинга.
    """
    p, d, loss_factor, angle_offset = opt_params
    theo_linear = compute_theoretical_grid_2d(angles, freqs, p, d, loss_factor, angle_offset, t_noise)
    theo_db = 10 * np.log10(np.maximum(theo_linear, 1e-12))
    
    residuals_linear = (exp_linear - theo_linear) * 100
    residuals_db = exp_db - theo_db
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.25)
    
    F, A = np.meshgrid(freqs, angles)
    
    # 1. Эксперимент
    im1 = axs[0, 0].pcolormesh(F, A, exp_db, cmap='viridis', shading='auto', vmin=-40, vmax=0)
    fig.colorbar(im1, ax=axs[0, 0], label='Пропускание (дБ)')
    axs[0, 0].set_title('Эксперимент T(f, theta) [дБ]')
    axs[0, 0].set_ylabel('Угол (град)')
    
    # 2. Прецизионная модель
    im2 = axs[0, 1].pcolormesh(F, A, theo_db, cmap='viridis', shading='auto', vmin=-40, vmax=0)
    fig.colorbar(im2, ax=axs[0, 1], label='Пропускание (дБ)')
    axs[0, 1].set_title(f'Модель Blanco + Noise Floor [дБ]\n(P={p*1e6:.2f} мкм, D={d*1e6:.2f} мкм, Offset={angle_offset:.3f}°)')
    
    # 3. Линейные остатки
    v_max_lin = max(abs(np.min(residuals_linear)), abs(np.max(residuals_linear)), 1.0)
    im3 = axs[1, 0].pcolormesh(F, A, residuals_linear, cmap='RdBu_r', shading='auto', vmin=-v_max_lin, vmax=v_max_lin)
    fig.colorbar(im3, ax=axs[1, 0], label='Отклонение (%)')
    axs[1, 0].set_title('Линейные остатки (Эксперимент - Модель) [%]')
    axs[1, 0].set_xlabel('Частота (ТГц)')
    axs[1, 0].set_ylabel('Угол (град)')
    
    # 4. Логарифмические остатки
    v_max_db = max(abs(np.min(residuals_db)), abs(np.max(residuals_db)), 1.0)
    im4 = axs[1, 1].pcolormesh(F, A, residuals_db, cmap='RdBu_r', shading='auto', vmin=-v_max_db, vmax=v_max_db)
    fig.colorbar(im4, ax=axs[1, 1], label='Отклонение (дБ)')
    axs[1, 1].set_title('Логарифмические остатки [дБ]')
    axs[1, 1].set_xlabel('Частота (ТГц)')
    
    fig.suptitle(f"Прецизионный 2D-анализ и аппроксимация {title_suffix}", fontsize=16, y=0.98)
    filename = f"plots_2d_fit_{title_suffix.lower().replace(' ', '_')}.png"
    plt.savefig(Path(__file__).parent / filename, dpi=150, bbox_inches='tight')
    plt.close()

def plot_slice_comparison(angles, freqs, exp_linear, exp_db, t_noise, opt_params, title_suffix=""):
    """
    Строит графики срезов на фиксированных частотах (0.4, 0.8, 1.2 ТГц)
    для детального сравнения теории и эксперимента.
    """
    p, d, loss_factor, angle_offset = opt_params
    theo_linear = compute_theoretical_grid_2d(angles, freqs, p, d, loss_factor, angle_offset, t_noise)
    theo_db = 10 * np.log10(np.maximum(theo_linear, 1e-12))
    
    target_freqs = [0.4, 0.8, 1.2]
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, f_target in enumerate(target_freqs):
        f_idx = np.argmin(np.abs(freqs - f_target))
        actual_f = freqs[f_idx]
        
        axs[idx].plot(angles, exp_db[:, f_idx], 'ko', label='Эксперимент')
        axs[idx].plot(angles, theo_db[:, f_idx], 'r-', linewidth=2, label='Прецизионная модель')
        
        # Нарисуем также паспортную идеальную модель для сравнения
        pas_theo_lin = compute_theoretical_grid_2d(angles, freqs, config.P_DEFAULT, config.D_DEFAULT, 0.15, 0.0, np.zeros_like(t_noise))
        pas_theo_db = 10 * np.log10(np.maximum(pas_theo_lin, 1e-12))
        axs[idx].plot(angles, pas_theo_db[:, f_idx], 'b:', label='Паспортная Blanco')
        
        axs[idx].set_title(f'Срез на частоте {actual_f:.2f} ТГц')
        axs[idx].set_xlabel('Угол (град)')
        axs[idx].set_ylabel('Пропускание (дБ)')
        axs[idx].set_ylim(-45, 2)
        axs[idx].grid(True, linestyle=':', alpha=0.6)
        axs[idx].legend(loc='lower left')
        
    fig.suptitle(f"Срезы угловой зависимости {title_suffix}", fontsize=16, y=0.98)
    filename = f"plots_2d_slices_{title_suffix.lower().replace(' ', '_')}.png"
    plt.savefig(Path(__file__).parent / filename, dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    print("=== Прецизионная 2D-аппроксимация ТГц-спектров по Blanco (Серия 2) ===")
    
    # 1. Загрузка данных
    try:
        angles, freqs, exp_linear, exp_db, t_noise = load_2d_experimental_data(repetition=2)
        print(f"Успешно загружены экспериментальные данные Серии 2:")
        print(f"  - Количество углов: {len(angles)} ({min(angles)}° - {max(angles)}°)")
        print(f"  - Частотный диапазон: {min(freqs):.2f} ТГц - {max(freqs):.2f} ТГц")
        print(f"  - Спектральный шумовой пол: {10*np.log10(np.mean(t_noise)):.1f} дБ (средний)")
    except Exception as e:
        print(f"Ошибка при загрузке: {e}")
        sys.exit(1)
        
    # 2. Запуск оптимизации
    print("\nЗапуск Nelder-Mead совместного подбора параметров...")
    opt_params, info = run_2d_optimization(angles, freqs, exp_linear, exp_db, t_noise, w_linear=10.0, w_db=1.0)
    
    print(f"\nРезультаты оптимизации:")
    print(f"  Статус успеха:      {info['success']} ({info['message']})")
    print(f"  Оптимальный шаг P:  {info['p_um']:.3f} мкм (паспортный: 16.0 мкм)")
    print(f"  Оптимальный диам D: {info['d_um']:.3f} мкм (паспортный: 11.0 мкм)")
    print(f"  Коэффициент потерь: {info['loss_factor']:.4f} дБ/ТГц")
    print(f"  Сдвиг нуля углов:   {info['angle_offset']:.3f}°")
    print(f"  Время расчета:      {info['time_sec']:.2f} сек")
    print(f"  Итоговая RMSE Lin:  {info['rmse_lin_pct']:.3f}%")
    print(f"  Итоговая RMSE dB:   {info['rmse_db']:.3f} дБ")
    
    # 3. Визуализация
    plot_2d_results(angles, freqs, exp_linear, exp_db, t_noise, opt_params, title_suffix="Precision Fit")
    plot_slice_comparison(angles, freqs, exp_linear, exp_db, t_noise, opt_params, title_suffix="Precision Fit")
    
    print("\nГрафики успешно построены и сохранены в src/!")
