import numpy as np
import lmfit
import time
import logging
from unified_optimizer import config
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.utils import find_auto_water_mask

def residual_2d_complex(params, angles_val, analysis_freqs, exp_trans_2d, valid_mask):
    """
    Функция невязки для lmfit. Возвращает 1D массив взвешенных отклонений.
    """
    p_um = params['P_um'].value
    d_um = params['D_um'].value
    loss_factor = params['loss_factor'].value
    gamma = params['gamma'].value if 'gamma' in params else config.GAMMA_DEFAULT
    angle_offset = params['angle_offset'].value
    tau_ps = params['tau_ps'].value

    p = p_um * 1e-6
    d = d_um * 1e-6

    if d >= p:
        return np.ones(np.sum(valid_mask) * 2) * 1e6

    theo_complex = compute_theoretical_grid_2d(
        angles_val, analysis_freqs, p, d, loss_factor, angle_offset, tau_ps, gamma=gamma
    )

    exp_masked = exp_trans_2d[valid_mask]
    theo_masked = theo_complex[valid_mask]

    # Амплитудная невязка
    amp_residual = np.abs(exp_masked) - np.abs(theo_masked)
    
    # Фазовая невязка
    phase_residual = np.angle(exp_masked) - np.angle(theo_masked)
    phase_residual = np.arctan2(np.sin(phase_residual), np.cos(phase_residual))

    # Веса
    W_AMP = 1.0
    W_PHASE = np.sqrt(0.1)

    return np.concatenate([amp_residual * W_AMP, phase_residual * W_PHASE])

def run_lmfit_2d(data_dict, dataset_name=""):
    """
    Подготовка данных и запуск LMFIT. Возвращает объект результата (MinimizerResult).
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
            
    f_max = getattr(config, 'F_MAX_2D', config.F_MAX)
    target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[target_indices]
    
    exp_trans_2d = np.zeros((len(angles_val), len(analysis_freqs)), dtype=complex)
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]
        
    # Маска линий воды
    Y_mean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12)), axis=0)
    auto_water_mask, _ = find_auto_water_mask(analysis_freqs, Y_mean_db)
    
    valid_mask = np.ones_like(exp_trans_2d, dtype=bool)
    for i in range(len(angles_val)):
        valid_mask[i, ~auto_water_mask] = False
        
    # Инициализация параметров LMFIT
    params = lmfit.Parameters()
    p_fixed_um = (config.P_FIXED * 1e6) if config.P_FIXED is not None else (config.P_DEFAULT * 1e6)
    
    # Снимаем жесткие границы (оставляем только очень широкие для физической осмысленности)
    params.add('P_um', value=p_fixed_um, min=1.0, max=100.0, vary=(config.P_FIXED is None))
    params.add('D_um', value=config.D_DEFAULT * 1e6, min=0.1, max=50.0)
    
    init_loss = 0.3 if config.USE_POWER_LAW else 0.15
    params.add('loss_factor', value=init_loss, min=0.0, max=5.0)
    
    if config.USE_POWER_LAW:
        params.add('gamma', value=config.GAMMA_DEFAULT, min=0.1, max=3.0, vary=config.OPTIMIZE_GAMMA)
        
    params.add('angle_offset', value=0.0, min=-10.0, max=10.0)
    params.add('tau_ps', value=0.0, min=-10.0, max=10.0)
    
    logging.info(f"[{dataset_name}] Starting LMFIT optimization...")
    start_time = time.time()
    
    mini = lmfit.Minimizer(residual_2d_complex, params, fcn_args=(angles_val, analysis_freqs, exp_trans_2d, valid_mask))
    result = mini.minimize(method='leastsq')
    
    elapsed = time.time() - start_time
    logging.info(f"[{dataset_name}] LMFIT finished in {elapsed:.2f}s. success={result.success}, nfev={result.nfev}")
    
    return result, mini
