import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path
from scipy.optimize import minimize

# Add src path
src_dir = Path(r"C:\Users\pop\.\.gemini\antigravity\worktrees\THz-Spectroscopy-Python-Manual2\apply-latest-updates\src")
sys.path.append(str(src_dir))

import theoretical
import fitting_2d

def run_2d_power_law_fit_fixed_geom():
    print("=== 2D OPTIMIZATION FOR POWER-LAW LOSS MODEL (FIXED GEOMETRY) ===")
    
    # 1. Load experimental data (repetition 2)
    angles, freqs, T_exp_lin, T_exp_db, t_noise = fitting_2d.load_2d_experimental_data(
        repetition=2, freq_start=0.2, freq_end=1.5
    )
    
    c_light = 3e8
    
    # Fixed geometry from passport
    p = 15.50e-6
    d = 5.67e-6
    d_over_p = d / p
    
    # Compute base t_perp and t_par (independent of losses/offsets)
    t_perp_arr = []
    t_par_arr = []
    for f in freqs:
        lambda_m = c_light / (f * 1e12)
        p_over_lambda = p / lambda_m
        t_perp_arr.append(theoretical.compute_t_perp(p_over_lambda, d_over_p, N=15))
        t_par_arr.append(theoretical.compute_t_par(p_over_lambda, d_over_p, N=15))
        
    t_perp_arr = np.array(t_perp_arr)
    t_par_arr = np.array(t_par_arr)
    
    # 2. Define loss function for Nelder-Mead
    # params: [loss_factor_db, angle_offset_deg, loss_exponent]
    def loss_2d(params):
        loss_factor_db, angle_offset, loss_exponent = params
        
        # Physical bounds
        if loss_factor_db < 0.0 or loss_factor_db > 5.0:
            return 1e10
        if angle_offset < -5.0 or angle_offset > 5.0:
            return 1e10
        if loss_exponent < 0.1 or loss_exponent > 3.0:
            return 1e10
            
        loss_factor_np = loss_factor_db / 4.343
        
        adjusted_angles_rad = np.deg2rad(angles - angle_offset)
        cos_a = np.cos(adjusted_angles_rad)[:, np.newaxis]
        sin_a = np.sin(adjusted_angles_rad)[:, np.newaxis]
        
        t_perp = t_perp_arr[np.newaxis, :]
        t_par = t_par_arr[np.newaxis, :]
        
        # Apply power-law loss model
        t_perp_eff = t_perp * np.exp(-0.5 * loss_factor_np * (freqs ** loss_exponent))[np.newaxis, :]
        t_par_eff = t_par * np.exp(-0.5 * loss_factor_np * (freqs ** loss_exponent))[np.newaxis, :]
        
        E_out = cos_a**2 * t_perp_eff + sin_a**2 * t_par_eff
        T_ideal = np.abs(E_out)**2
        T_mod_lin = T_ideal + t_noise[np.newaxis, :]
        T_mod_lin = np.clip(T_mod_lin, 0.0, 20.0)
        T_mod_db = 10 * np.log10(np.maximum(T_mod_lin, 1e-12))
        
        # We calculate RMSE on valid points (above noise floor)
        valid_mask = T_exp_lin > 1.5 * t_noise[np.newaxis, :]
        return np.mean((T_exp_db[valid_mask] - T_mod_db[valid_mask]) ** 2)
        
    # 3. Initial guess
    # loss_factor_db=0.368, angle_offset=-0.45, loss_exponent=1.58
    init_guess = [0.368, -0.45, 1.58]
    
    print("Initial RMSE on valid points:", np.sqrt(loss_2d(init_guess)))
    
    res = minimize(loss_2d, init_guess, method='Nelder-Mead', options={'maxiter': 500})
    loss_factor_db_opt, offset_opt, exponent_opt = res.x
    rmse_opt = np.sqrt(res.fun)
    
    print("\n--- OPTIMIZATION RESULTS (FIXED GEOMETRY) ---")
    print(f"Fixed period P: {p*1e6:.4f} um")
    print(f"Fixed strip width D: {d*1e6:.4f} um")
    print(f"Optimal loss factor (dB): {loss_factor_db_opt:.4f} dB/THz^{exponent_opt:.2f}")
    print(f"Optimal loss exponent (gamma): {exponent_opt:.4f}")
    print(f"Optimal angular offset: {offset_opt:.4f} deg")
    print(f"Minimal 2D RMSE on valid points: {rmse_opt:.4f} dB")
    
    # 4. Generate residuals with optimal parameters
    loss_factor_np = loss_factor_db_opt / 4.343
    
    adjusted_angles_rad = np.deg2rad(angles - offset_opt)
    cos_a = np.cos(adjusted_angles_rad)[:, np.newaxis]
    sin_a = np.sin(adjusted_angles_rad)[:, np.newaxis]
    
    t_perp = t_perp_arr[np.newaxis, :]
    t_par = t_par_arr[np.newaxis, :]
    
    # Apply optimal power-law loss model
    t_perp_eff = t_perp * np.exp(-0.5 * loss_factor_np * (freqs ** exponent_opt))[np.newaxis, :]
    t_par_eff = t_par * np.exp(-0.5 * loss_factor_np * (freqs ** exponent_opt))[np.newaxis, :]
    
    E_out = cos_a**2 * t_perp_eff + sin_a**2 * t_par_eff
    T_ideal = np.abs(E_out)**2
    T_mod_lin = T_ideal + t_noise[np.newaxis, :]
    T_mod_lin = np.clip(T_mod_lin, 0.0, 20.0)
    T_mod_db = 10 * np.log10(np.maximum(T_mod_lin, 1e-12))
    
    residuals_db = T_exp_db - T_mod_db
    
    # 5. Plot heatmap of residuals
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(residuals_db, aspect='auto', extent=[freqs[0], freqs[-1], angles[-1], angles[0]], 
               cmap='RdBu_r', vmin=-3.0, vmax=3.0)
    fig.colorbar(im, label='Residual (Experiment - Model) [dB]')
    ax.invert_yaxis()
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Rotator Angle (deg)')
    ax.set_title(f'2D Residuals Heatmap (RMSE: {rmse_opt:.3f} dB)\nFixed Geometry: P={p*1e6:.2f} um, D={d*1e6:.2f} um')
    
    output_path = Path(r"C:\Users\pop\.gemini\antigravity\brain\23ae18f8-391f-4212-8240-c4d5640962ff\residuals_heatmap_2d_fixed_geom.png")
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved fixed-geometry residuals heatmap to: {output_path}")

if __name__ == '__main__':
    run_2d_power_law_fit_fixed_geom()
