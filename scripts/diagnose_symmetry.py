import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from unified_optimizer.data_manager import DataManager
from unified_optimizer import config

def diagnose():
    manager = DataManager(config.DATA_DIR)
    for ds in ['series1', 'series2']:
        print(f"\n--- Диагностика {ds} ---")
        try:
            data = manager.get_data_for_dataset(ds)
        except Exception:
            continue
            
        angles = sorted(list(data.keys()))
        positive_angles = [a for a in angles if 0 <= a <= 90]
        
        for pos_a in positive_angles:
            neg_a = -pos_a
            if neg_a in data and pos_a != 0:
                t_s_pos, E_s_pos, _, _ = data[pos_a]
                t_s_neg, E_s_neg, _, _ = data[neg_a]
                
                # Сравниваем интегралы энергии
                int_pos = np.trapezoid(E_s_pos**2, t_s_pos)
                int_neg = np.trapezoid(E_s_neg**2, t_s_neg)
                
                diff_percent = abs(int_pos - int_neg) / ((int_pos + int_neg)/2) * 100
                print(f"Угол {pos_a:5.1f} vs {neg_a:5.1f}: разница {diff_percent:.1f}%")

if __name__ == "__main__":
    diagnose()
