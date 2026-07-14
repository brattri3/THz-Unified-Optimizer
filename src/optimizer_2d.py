import numpy as np
from scipy.optimize import minimize
import time

from unified_optimizer import config, model_blanco

def get_transmission_spectra(t_s, E_s, t_b, E_b):
    """
    Рассчитывает спектр пропускания и спектр фона для 2D-анализа.
    """
    dt = t_s[1] - t_s[0]
    N = len(t_s)
    fft_s = np.fft.rfft(E_s)
    fft_b = np.fft.rfft(E_b)
    freq = np.fft.rfftfreq(N, dt)
    
    spec_s = np.abs(fft_s)
    spec_b = np.abs(fft_b)
    
    # Избегаем деления на ноль
    transmission = np.zeros_like(spec_b)
    valid = spec_b > 1e-10
    transmission[valid] = (spec_s[valid] / spec_b[valid])**2
    
    return freq, spec_s, spec_b, transmission

def compute_theoretical_grid_2d(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, t_noise, N=15):
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
    
    adjusted_angles_rad = np.deg2rad(angles_deg - angle_offset)
    cos_a = np.cos(adjusted_angles_rad)[:, np.newaxis]
    sin_a = np.sin(adjusted_angles_rad)[:, np.newaxis]
    
    t_perp = t_perp_arr[np.newaxis, :]
    t_par = t_par_arr[np.newaxis, :]
    
    t_perp_eff = t_perp * np.exp(-0.5 * loss_factor * freqs_thz)[np.newaxis, :]
    t_par_eff = t_par * np.exp(-0.5 * loss_factor * freqs_thz)[np.newaxis, :]
    
    E_out = cos_a**2 * t_perp_eff + sin_a**2 * t_par_eff
    T_ideal = np.abs(E_out)**2
    T_mod = T_ideal + t_noise[np.newaxis, :]
    
    return np.clip(T_mod, 0.0, 20.0)

def optimize_2d_spectral(data_dict):
    """
    Выполняет 2D спектрально-угловую оптимизацию (по Nelder-Mead).
    """
    angles = sorted(list(data_dict.keys()))
    angles_val = np.array(angles)
    
    bg_spectra = []
    individual_results = {}
    freqs_common = None
    
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        bg_spectra.append(s_b)
        individual_results[a] = trans
        if freqs_common is None:
            freqs_common = f
            
    bg_avg = np.mean(bg_spectra, axis=0)
    tail_indices = np.where(freqs_common >= 3.0)[0]
    if len(tail_indices) == 0:
        tail_indices = np.where(freqs_common >= 2.5)[0]
    noise_floor_amplitude = np.mean(bg_avg[tail_indices])
    dr_power = (bg_avg / noise_floor_amplitude)**2
    t_noise_floor = 1.0 / dr_power
    
    target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= config.F_MAX))[0]
    analysis_freqs = freqs_common[target_indices]
    fit_t_noise = t_noise_floor[target_indices]
    
    exp_trans_2d = np.zeros((len(angles_val), len(analysis_freqs)))
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]
        
    exp_trans_db_2d = 10 * np.log10(np.maximum(exp_trans_2d, 1e-12))
    
    def loss_function_2d(params_scaled):
        p_um, d_um, loss_factor, angle_offset = params_scaled
        
        if p_um <= 5.0 or p_um >= 40.0: return 1e8
        if d_um <= 1.0 or d_um >= 25.0: return 1e8
        if d_um >= p_um: return 1e8 + (d_um - p_um) * 1e9
        if loss_factor < 0.0 or loss_factor > 5.0: return 1e8
        if angle_offset < -10.0 or angle_offset > 10.0: return 1e8
            
        p = p_um * 1e-6
        d = d_um * 1e-6
            
        theo_linear = compute_theoretical_grid_2d(angles_val, analysis_freqs, p, d, loss_factor, angle_offset, fit_t_noise)
        theo_db = 10 * np.log10(np.maximum(theo_linear, 1e-12))
        
        valid_mask = exp_trans_2d > 1.5 * fit_t_noise[np.newaxis, :]
        if np.sum(valid_mask) == 0: return 1e8
            
        diff_lin = exp_trans_2d[valid_mask] - theo_linear[valid_mask]
        rmse_lin = np.sqrt(np.mean(diff_lin**2))
        rmse_db = np.sqrt(np.mean((exp_trans_db_2d[valid_mask] - theo_db[valid_mask])**2))
        
        return config.W_LINEAR * rmse_lin + config.W_DB * rmse_db

    initial_guess = [config.P_DEFAULT * 1e6, config.D_DEFAULT * 1e6, 0.0, 0.0]
    res = minimize(
        loss_function_2d, 
        initial_guess, 
        method='Nelder-Mead',
        options={'maxiter': 50, 'xatol': 1e-2, 'fatol': 1e-2}
    )
    
    p_f, d_f, a_f, th_f = res.x
    return {
        'P_eff_um': p_f,
        'D_eff_um': d_f,
        'loss_factor': a_f,
        'theta_offset': th_f,
        'success': res.success,
        'fun': res.fun
    }
