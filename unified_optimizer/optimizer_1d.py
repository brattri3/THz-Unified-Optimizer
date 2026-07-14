import numpy as np
from scipy.optimize import minimize
import time

from unified_optimizer import config, model_blanco

def build_spectral_basis(t_bg: np.ndarray, E_bg: np.ndarray):
    """Строит массив частот и весовой спектр |E_bg(ν)|² из фонового сигнала."""
    dt = t_bg[1] - t_bg[0]
    N = len(t_bg)
    fft = np.fft.rfft(E_bg)
    freq = np.fft.rfftfreq(N, dt)
    mask = (freq >= config.F_MIN) & (freq <= config.F_MAX)
    freqs = freq[mask]
    weights = np.abs(fft[mask])**2
    return freqs, weights

def theory_T_integral(angle_deg, p_eff_m, d_eff_m, alpha, gamma, freqs_THz, weights, theta_offset_deg=0.0, eps_floor=0.0):
    theta_rad = np.deg2rad(angle_deg + theta_offset_deg)
    T_blanco = np.empty(len(freqs_THz))
    for i, nu in enumerate(freqs_THz):
        lam = model_blanco.C_LIGHT / (nu * 1e12)
        p_over_lam = p_eff_m / lam
        d_over_p   = d_eff_m / p_eff_m
        T_blanco[i] = model_blanco.transmission_two_polarizers(theta_rad, p_over_lam, d_over_p, N=10)
    loss = np.exp(-alpha * freqs_THz**gamma)
    T_model = T_blanco * loss
    numerator = np.trapezoid(T_model * weights, freqs_THz)
    denominator = np.trapezoid(weights, freqs_THz)
    return float(numerator / denominator) + eps_floor

def build_theory_curve(angles, p_eff_m, d_eff_m, alpha, gamma, freqs_THz, weights, theta_offset_deg=0.0, eps_floor=0.0):
    return np.array([
        theory_T_integral(a, p_eff_m, d_eff_m, alpha, gamma, freqs_THz, weights, theta_offset_deg, eps_floor)
        for a in angles
    ])

def optimize_1d_integral(data_dict):
    """
    Выполняет 3-этапную 1D интегральную оптимизацию.
    Возвращает словарь с найденными параметрами.
    """
    angles = np.array(sorted(data_dict.keys()))
    T_exp = []
    bg_spectra_weights = []
    freqs_common = None
    
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        E_s_int = np.trapezoid(E_s**2, t_s)
        E_b_int = np.trapezoid(E_b**2, t_b)
        T_exp.append(E_s_int / E_b_int)
        
        f, w = build_spectral_basis(t_b, E_b)
        bg_spectra_weights.append(w)
        if freqs_common is None:
            freqs_common = f
            
    T_exp = np.array(T_exp)
    weights_avg = np.mean(bg_spectra_weights, axis=0)
    
    # ─── ЭТАП 1: Линейная шкала (4 параметра)
    def loss_stage1(x):
        p_um, d_um, alpha, gamma = x
        T_mod = build_theory_curve(angles, p_um*1e-6, d_um*1e-6, alpha, gamma, freqs_common, weights_avg)
        return np.sum((T_mod - T_exp)**2)
        
    x0_1 = [config.P_DEFAULT*1e6, config.D_DEFAULT*1e6, 0.1, 1.0]
    b1 = [(15.5, 16.5), (1.0, 15.0), (0.0, 5.0), (0.1, 3.0)]
    res1 = minimize(loss_stage1, x0_1, bounds=b1, method='L-BFGS-B')
    
    # ─── ЭТАП 2: дБ-шкала (4 параметра)
    def loss_stage2(x):
        p_um, d_um, alpha, gamma = x
        T_mod = build_theory_curve(angles, p_um*1e-6, d_um*1e-6, alpha, gamma, freqs_common, weights_avg)
        diff = 10 * np.log10(np.maximum(T_exp, 1e-12)) - 10 * np.log10(np.maximum(T_mod, 1e-12))
        return np.sum(diff**2)
        
    p_opt, d_opt, a_opt, g_opt = res1.x
    x0_2 = [p_opt, d_opt, a_opt, g_opt]
    b2 = [(p_opt-0.5, p_opt+0.5), (d_opt-0.5, d_opt+0.5), (0.0, 5.0), (0.1, 3.0)]
    res2 = minimize(loss_stage2, x0_2, bounds=b2, method='L-BFGS-B')
    
    # ─── ЭТАП 3: дБ + аппаратные параметры (6 параметров)
    def loss_stage3(x):
        p_um, d_um, alpha, gamma, th_off, eps = x
        T_mod = build_theory_curve(angles, p_um*1e-6, d_um*1e-6, alpha, gamma, freqs_common, weights_avg, th_off, eps)
        diff = 10 * np.log10(np.maximum(T_exp, 1e-12)) - 10 * np.log10(np.maximum(T_mod, 1e-12))
        return np.sum(diff**2)
        
    p_opt, d_opt, a_opt, g_opt = res2.x
    x0_3 = [p_opt, d_opt, a_opt, g_opt, 0.0, 0.0]
    b3 = [(p_opt-0.2, p_opt+0.2), (d_opt-0.2, d_opt+0.2), (0.0, 5.0), (0.1, 3.0), (-5.0, 5.0), (0.0, 1e-2)]
    res3 = minimize(loss_stage3, x0_3, bounds=b3, method='L-BFGS-B')
    
    p_f, d_f, a_f, g_f, th_f, eps_f = res3.x
    return {
        'P_eff_um': p_f,
        'D_eff_um': d_f,
        'alpha': a_f,
        'gamma': g_f,
        'theta_offset': th_f,
        'eps_floor': eps_f,
        'success': res3.success,
        'fun': res3.fun
    }
