import sys
import numpy as np
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from unified_optimizer.data_manager import DataManager
from unified_optimizer import config
from unified_optimizer.optimizer_1d import theory_T_integral, build_spectral_basis

P_OPT = 15.50e-6
D_OPT = 4.045e-6
OFFSET_OPT = 0.35
LOSS_OPT = 0.295
GAMMA_OPT = 1.69

def main():
    manager = DataManager(Path(config.DATA_DIR))
    datasets = manager.get_datasets()
    
    alpha_scaled = 0.5 * (LOSS_OPT / 4.343)
    
    for ds in datasets:
        data_dict = manager.get_data_for_dataset(ds)
        angles = np.array(sorted(list(data_dict.keys())))
        if len(angles) == 0: continue
        idx = np.argmin(np.abs(angles))
        a0 = angles[idx]
        
        t_s, E_s, t_b, E_b = data_dict[a0]
        
        E_s_clean = E_s - np.mean(E_s)
        E_b_clean = E_b - np.mean(E_b)
        
        E_s_int = np.trapezoid(E_s_clean**2, t_s)
        E_b_int = np.trapezoid(E_b_clean**2, t_b)
        T_exp = np.sqrt(E_s_int / E_b_int)
        
        f, w = build_spectral_basis(t_b, E_b_clean)
        T_theo = theory_T_integral(a0, P_OPT, D_OPT, alpha_scaled, GAMMA_OPT, f, w, OFFSET_OPT)
        
        diff_pct = (T_exp - T_theo) / T_theo * 100
        print(f"{ds:<10} | Exp: {T_exp:.4f} | Theo (Scaled): {T_theo:.4f} | Diff: {diff_pct:.2f}%")

if __name__ == "__main__":
    main()
