import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.stats as stats

sys.path.append(str(Path(__file__).resolve().parent.parent))
from unified_optimizer.data_manager import DataManager
from unified_optimizer import config
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.utils import find_auto_water_mask

# Parameters
P_OPT = 17.50e-6
D_OPT = 5.053e-6
OFFSET_OPT = -0.183
LOSS_OPT = 0.287
GAMMA_OPT = 1.868
TAU_PS_OPT = 0.0327

def get_residuals_for_ds(manager, ds):
    data_dict = manager.get_data_for_dataset(ds)
    angles = sorted(list(data_dict.keys()))
    
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
            
    bg_avg = np.mean(np.abs(bg_spectra), axis=0)
    tail_indices = np.where(freqs_common >= 3.0)[0]
    if len(tail_indices) == 0:
        tail_indices = np.where(freqs_common >= 2.5)[0]
    noise_floor_amplitude = np.mean(bg_avg[tail_indices])
    
    target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= config.F_MAX))[0]
    analysis_freqs = freqs_common[target_indices]
    
    exp_trans_2d = np.zeros((len(angles), len(analysis_freqs)), dtype=np.complex128)
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]
        
    Y_mean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12)), axis=0)
    water_mask, _ = find_auto_water_mask(analysis_freqs, Y_mean_db)
        
    theo_trans_2d = compute_theoretical_grid_2d(
        angles, analysis_freqs, P_OPT, D_OPT, 
        LOSS_OPT, OFFSET_OPT, TAU_PS_OPT, gamma=GAMMA_OPT
    )
    
    valid_2d_mask = np.abs(exp_trans_2d) > 2.0 * (noise_floor_amplitude / np.abs(np.mean(bg_spectra, axis=0)[target_indices]))
    valid_2d_mask = valid_2d_mask & water_mask[np.newaxis, :]
    
    amp_res = 20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12)) - 20 * np.log10(np.maximum(np.abs(theo_trans_2d), 1e-12))
    phase_res = np.angle(exp_trans_2d) - np.angle(theo_trans_2d)
    phase_res = (phase_res + np.pi) % (2 * np.pi) - np.pi
    
    return angles, analysis_freqs, valid_2d_mask, amp_res, phase_res

def main():
    manager = DataManager(config.DATA_DIR)
    
    ang1, f1, mask1, amp1, ph1 = get_residuals_for_ds(manager, '356att')
    ang2, f2, mask2, amp2, ph2 = get_residuals_for_ds(manager, 'series3')
    
    common_angles = sorted(list(set(ang1).intersection(set(ang2))))
    
    if not np.allclose(f1, f2):
        print("Frequencies do not match!")
        return
        
    ang_list = []
    corr_a_list = []
    corr_p_list = []
        
    for ang in common_angles:
        idx1 = ang1.index(ang)
        idx2 = ang2.index(ang)
        
        m_common = mask1[idx1] & mask2[idx2]
        a1 = amp1[idx1][m_common]
        a2 = amp2[idx2][m_common]
        p1 = ph1[idx1][m_common]
        p2 = ph2[idx2][m_common]
        
        if len(a1) > 10:
            corr_a, _ = stats.pearsonr(a1, a2)
            corr_p, _ = stats.pearsonr(p1, p2)
            ang_list.append(ang)
            corr_a_list.append(corr_a)
            corr_p_list.append(corr_p)

    a1_all = []
    a2_all = []
    p1_all = []
    p2_all = []
    
    for ang in common_angles:
        idx1 = ang1.index(ang)
        idx2 = ang2.index(ang)
        m_common = mask1[idx1] & mask2[idx2]
        a1_all.extend(amp1[idx1][m_common])
        a2_all.extend(amp2[idx2][m_common])
        p1_all.extend(ph1[idx1][m_common])
        p2_all.extend(ph2[idx2][m_common])

    # Plot 1: Angular dependence of correlation
    plt.figure(figsize=(10, 6))
    plt.plot(ang_list, corr_a_list, 'o-', label='Amplitude Correlation $R_A$', color='red')
    plt.plot(ang_list, corr_p_list, 's-', label='Phase Correlation $R_\phi$', color='blue')
    plt.axhline(0, color='gray', linestyle='--')
    plt.xlabel('Angle (deg)')
    plt.ylabel('Pearson Correlation Coefficient (R)')
    plt.title('Angular Dependence of Residual Correlations (356att vs series3)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(Path(config.BASE_DIR).parent / "docs" / "images" / "correlation_angular.png")
    plt.close()

    # Plot 2: Scatter plot of global residuals
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.scatter(a1_all, a2_all, alpha=0.3, color='red', s=5)
    plt.xlabel('Amplitude Residuals (356att) [dB]')
    plt.ylabel('Amplitude Residuals (series3) [dB]')
    plt.title(f'Amplitude Global Corr (R={stats.pearsonr(a1_all, a2_all)[0]:.2f})')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.scatter(p1_all, p2_all, alpha=0.3, color='blue', s=5)
    plt.xlabel('Phase Residuals (356att) [rad]')
    plt.ylabel('Phase Residuals (series3) [rad]')
    plt.title(f'Phase Global Corr (R={stats.pearsonr(p1_all, p2_all)[0]:.2f})')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(Path(config.BASE_DIR).parent / "docs" / "images" / "correlation_scatter.png")
    plt.close()

    print(f"Global Amplitude Correlation R: {stats.pearsonr(a1_all, a2_all)[0]:.4f}")
    print(f"Global Phase Correlation R: {stats.pearsonr(p1_all, p2_all)[0]:.4f}")
    print("Plots saved to docs/images/")
    
if __name__ == "__main__":
    main()
