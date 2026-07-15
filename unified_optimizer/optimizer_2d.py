import numpy as np
from scipy.optimize import minimize
import time

from unified_optimizer import config, model_blanco

def get_transmission_spectra(t_s, E_s, t_b, E_b):
    """
    Рассчитывает спектр пропускания и спектр фона для 2D-анализа.
    Теперь возвращает комплексное отношение (амплитуду и фазу).
    """
    dt = t_s[1] - t_s[0]
    N = len(t_s)
    fft_s = np.fft.rfft(E_s)
    fft_b = np.fft.rfft(E_b)
    freq = np.fft.rfftfreq(N, dt)
    
    spec_s = np.abs(fft_s)
    spec_b = np.abs(fft_b)
    
    # Избегаем деления на ноль
    transmission = np.zeros_like(fft_b, dtype=complex)
    valid = spec_b > 1e-10
    transmission[valid] = fft_s[valid] / fft_b[valid]
    
    return freq, spec_s, spec_b, transmission

def compute_theoretical_grid_2d(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps, gamma=1.0, N=15):
    """
    Возвращает комплексный массив пропускания модели.
    """
    t_perp_arr = []
    t_par_arr = []
    d_over_p = d / p
    for f in freqs_thz:
        lambda_m = model_blanco.C_LIGHT / (f * 1e12)
        p_over_lambda = p / lambda_m
        t_perp_arr.append(model_blanco.compute_t_perp(p_over_lambda, d_over_p, N))
        t_par_arr.append(model_blanco.compute_t_par(p_over_lambda, d_over_p, N))
    
    t_perp_arr = np.array(t_perp_arr)
    t_par_arr = np.array(t_par_arr)
    
    angles_deg = np.array(angles_deg)
    adjusted_angles_rad = np.deg2rad(angles_deg - angle_offset)
    cos_a = np.cos(adjusted_angles_rad)[:, np.newaxis]
    sin_a = np.sin(adjusted_angles_rad)[:, np.newaxis]
    
    t_perp = t_perp_arr[np.newaxis, :]
    t_par = t_par_arr[np.newaxis, :]
    
    if config.USE_POWER_LAW:
        loss_factor_np = loss_factor / 4.343
        loss_amp = np.exp(-0.5 * loss_factor_np * (freqs_thz ** gamma))[np.newaxis, :]
    else:
        loss_amp = np.exp(-0.5 * loss_factor * freqs_thz)[np.newaxis, :]
        
    t_perp_eff = t_perp * loss_amp
    t_par_eff = t_par * loss_amp
    
    # Комплексная матрица Джонса
    E_out = cos_a**2 * t_perp_eff + sin_a**2 * t_par_eff
    
    # Фазовая задержка из-за разности оптического пути
    phase_delay = np.exp(-1j * 2 * np.pi * freqs_thz * tau_ps)[np.newaxis, :]
    
    return E_out * phase_delay

def optimize_2d_spectral(data_dict):
    """
    Выполняет прецизионную 2D спектрально-угловую оптимизацию (по Nelder-Mead).
    Использует комплексную целевую функцию и маскирует линии воды.
    """
    angles = sorted(list(data_dict.keys()))
    if config.ANGLES_LIMIT_2D is not None:
        min_a, max_a = config.ANGLES_LIMIT_2D
        angles = [a for a in angles if min_a <= a <= max_a]
        
    angles_val = np.array(angles)
    
    bg_spectra = []
    individual_results = {}
    freqs_common = None
    
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans_complex = get_transmission_spectra(t_s, E_s, t_b, E_b)
        bg_spectra.append(s_b)
        individual_results[a] = trans_complex
        if freqs_common is None:
            freqs_common = f
            
    bg_avg = np.mean(bg_spectra, axis=0)
    tail_indices = np.where(freqs_common >= 3.0)[0]
    if len(tail_indices) == 0:
        tail_indices = np.where(freqs_common >= 2.5)[0]
    noise_floor_amplitude = np.mean(bg_avg[tail_indices])
    dr_power = (bg_avg / noise_floor_amplitude)**2
    t_noise_floor = 1.0 / dr_power
    
    f_max = getattr(config, 'F_MAX_2D', config.F_MAX)
    target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[target_indices]
    fit_t_noise = t_noise_floor[target_indices]
    
    exp_trans_2d = np.zeros((len(angles_val), len(analysis_freqs)), dtype=complex)
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]
        
    free_params = []
    if config.P_FIXED is None:
        free_params.append(('P_um', config.P_DEFAULT * 1e6, (5.0, 40.0)))
    free_params.append(('D_um', config.D_DEFAULT * 1e6, (1.0, 25.0)))
    
    init_loss = 0.3 if config.USE_POWER_LAW else 0.15
    free_params.append(('loss_factor', init_loss, (0.0, 5.0)))
    
    if config.USE_POWER_LAW and config.OPTIMIZE_GAMMA:
        free_params.append(('gamma', config.GAMMA_DEFAULT, (0.1, 3.0)))
    free_params.append(('angle_offset', 0.0, (-10.0, 10.0)))
    free_params.append(('tau_ps', 0.0, (-10.0, 10.0))) # Фазовая задержка
    
    initial_guess = [p[1] for p in free_params]
    
    def loss_function_2d(free_vals):
        p_um = config.P_FIXED * 1e6 if config.P_FIXED is not None else 16.0
        d_um = config.D_DEFAULT * 1e6
        loss_factor = 0.3 if config.USE_POWER_LAW else 0.15
        gamma = config.GAMMA_DEFAULT if config.USE_POWER_LAW else 1.0
        angle_offset = 0.0
        tau_ps = 0.0
        
        for val, (name, _, bounds) in zip(free_vals, free_params):
            if val < bounds[0] or val > bounds[1]:
                return 1e8
            if name == 'P_um': p_um = val
            elif name == 'D_um': d_um = val
            elif name == 'loss_factor': loss_factor = val
            elif name == 'gamma': gamma = val
            elif name == 'angle_offset': angle_offset = val
            elif name == 'tau_ps': tau_ps = val
            
        if d_um >= p_um:
            return 1e8 + (d_um - p_um) * 1e9
            
        p = p_um * 1e-6
        d = d_um * 1e-6
        
        theo_complex = compute_theoretical_grid_2d(angles_val, analysis_freqs, p, d, loss_factor, angle_offset, tau_ps, gamma=gamma)
        
        # Маска достоверных точек
        valid_mask = np.abs(exp_trans_2d) > 1.5 * np.sqrt(fit_t_noise[np.newaxis, :])
        
        # Вырезание линий воды
        if hasattr(config, 'WATER_LINES_THZ'):
            for w_min, w_max in config.WATER_LINES_THZ:
                water_idx = (analysis_freqs >= w_min) & (analysis_freqs <= w_max)
                valid_mask[:, water_idx] = False
                
        if np.sum(valid_mask) == 0:
            return 1e6
            
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

    res = minimize(
        loss_function_2d, 
        initial_guess, 
        method='Nelder-Mead',
        options={'maxiter': 1000, 'xatol': 1e-4, 'fatol': 1e-4}
    )
    
    final_vals = res.x
    p_um = config.P_FIXED * 1e6 if config.P_FIXED is not None else 16.0
    d_um = config.D_DEFAULT * 1e6
    loss_factor = 0.3 if config.USE_POWER_LAW else 0.15
    gamma = config.GAMMA_DEFAULT if config.USE_POWER_LAW else 1.0
    angle_offset = 0.0
    tau_ps = 0.0
    
    for val, (name, _, _) in zip(final_vals, free_params):
        if name == 'P_um': p_um = val
        elif name == 'D_um': d_um = val
        elif name == 'loss_factor': loss_factor = val
        elif name == 'gamma': gamma = val
        elif name == 'angle_offset': angle_offset = val
        elif name == 'tau_ps': tau_ps = val
        
    return {
        'P_eff_um': p_um,
        'D_eff_um': d_um,
        'loss_factor': loss_factor,
        'gamma': gamma,
        'theta_offset': angle_offset,
        'tau_ps': tau_ps,
        'success': res.success,
        'fun': res.fun
    }
