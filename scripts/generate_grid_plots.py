import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from unified_optimizer.data_manager import DataManager
from unified_optimizer import config
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.optimizer_1d import theory_T_integral, build_spectral_basis

# Используем оптимальные параметры из 2D оптимизации
P_OPT = 15.50e-6
D_OPT = 4.045e-6
OFFSET_OPT = 0.35
LOSS_OPT = 0.295
GAMMA_OPT = 1.69
TAU_PS_OPT = 0.029

def main():
    data_dir = Path(config.DATA_DIR)
    images_dir = Path(config.BASE_DIR.parent / "docs" / "images")
    images_dir.mkdir(exist_ok=True, parents=True)
    
    manager = DataManager(data_dir)
    datasets = manager.get_datasets()
    
    report_path = Path(config.BASE_DIR.parent / "docs" / "artifacts" / "grid_plots_report.md")
    report_path.parent.mkdir(exist_ok=True, parents=True)
    
    with open(report_path, "w", encoding="utf-8") as f_report:
        f_report.write("# Угловые зависимости и невязки: Спектральный и Интегральный методы\n\n")
        f_report.write("## Описание модели и используемые параметры\n\n")
        f_report.write("В данном отчете представлено сопоставление экспериментальных данных ТГц-TDS спектроскопии с теоретической электродинамической моделью Бланко для системы из двух проволочных поляризаторов (аттенюатора).\n\n")
        f_report.write("### Оптимизированные параметры решётки и тракта:\n")
        f_report.write(f"- **Период проволочной решётки ($P$):** {P_OPT * 1e6:.3f} мкм\n")
        f_report.write(f"- **Диаметр проволоки ($D$):** {D_OPT * 1e6:.3f} мкм\n")
        f_report.write(f"- **Систематический сдвиг угла ($\\theta_{{offset}}$):** {OFFSET_OPT:.2f}^\\circ\n")
        f_report.write(f"- **Коэффициент поглощения/рассеяния ($loss\\_factor$):** {LOSS_OPT:.3f} (степень $\\gamma = {GAMMA_OPT}$)\n")
        f_report.write(f"- **Фазовая задержка оптического пути ($\\tau_{{ps}}$):** {TAU_PS_OPT:.3f} пс\n\n")
        f_report.write("### Описание методов:\n")
        f_report.write(r"1. **Спектральный анализ:** Амплитудный коэффициент пропускания $T(\nu) = |E_s(\nu)| / |E_b(\nu)|$ оценивается на фиксированных частотах $0.5$ ТГц и $1.0$ ТГц." + "\n")
        f_report.write(r"2. **Интегральный метод (Time-Domain):** Амплитудное пропускание оценивается по полной энергии импульса после удаления DC-составляющей (смещения ноля): $T_{int} = \sqrt{\int E_s^2(t) dt / \int E_b^2(t) dt}$. Теоретическое значение вычисляется как RMS-усреднение спектральной модели по фоновому весовому спектру." + "\n\n")
        f_report.write("---\n\n")
        
        for ds in datasets:
            print(f"Обработка {ds}...")
            data_dict = manager.get_data_for_dataset(ds)
            angles = np.array(sorted(list(data_dict.keys())))
            if len(angles) == 0: continue
            
            # --- 1. Spectral Analysis ---
            individual_results = {}
            freqs_common = None
            bg_spectra = []
            
            for a in angles:
                t_s, E_s, t_b, E_b = data_dict[a]
                f_axis, s_s, s_b, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
                individual_results[a] = trans
                bg_spectra.append(s_b)
                if freqs_common is None:
                    freqs_common = f_axis
                    
            target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= config.F_MAX))[0]
            analysis_freqs = freqs_common[target_indices]
            
            idx_05 = np.argmin(np.abs(analysis_freqs - 0.5))
            idx_10 = np.argmin(np.abs(analysis_freqs - 1.0))
            
            f_05 = analysis_freqs[idx_05]
            f_10 = analysis_freqs[idx_10]
            
            exp_trans_2d = np.zeros((len(angles), len(analysis_freqs)), dtype=np.complex128)
            for i, ang in enumerate(angles):
                exp_trans_2d[i, :] = individual_results[ang][target_indices]
                
            theo_trans_2d = compute_theoretical_grid_2d(
                angles, analysis_freqs, P_OPT, D_OPT, 
                LOSS_OPT, OFFSET_OPT, TAU_PS_OPT, gamma=GAMMA_OPT
            )
            
            # Создаем плотную сетку углов от -90 до 90 для гладкой теоретической кривой
            angles_dense = np.linspace(-90, 90, 361)
            
            theo_trans_2d_dense = compute_theoretical_grid_2d(
                angles_dense, analysis_freqs, P_OPT, D_OPT, 
                LOSS_OPT, OFFSET_OPT, TAU_PS_OPT, gamma=GAMMA_OPT
            )
            
            exp_amp_05 = np.abs(exp_trans_2d[:, idx_05])
            theo_amp_05_exp = np.abs(theo_trans_2d[:, idx_05])
            theo_amp_05_dense = np.abs(theo_trans_2d_dense[:, idx_05])
            
            exp_amp_10 = np.abs(exp_trans_2d[:, idx_10])
            theo_amp_10_exp = np.abs(theo_trans_2d[:, idx_10])
            theo_amp_10_dense = np.abs(theo_trans_2d_dense[:, idx_10])
            
            res_lin_05 = exp_amp_05 - theo_amp_05_exp
            res_lin_10 = exp_amp_10 - theo_amp_10_exp
            
            exp_db_05 = 20 * np.log10(np.maximum(exp_amp_05, 1e-12))
            theo_db_05_exp = 20 * np.log10(np.maximum(theo_amp_05_exp, 1e-12))
            theo_db_05_dense = 20 * np.log10(np.maximum(theo_amp_05_dense, 1e-12))
            res_db_05 = exp_db_05 - theo_db_05_exp
            
            exp_db_10 = 20 * np.log10(np.maximum(exp_amp_10, 1e-12))
            theo_db_10_exp = 20 * np.log10(np.maximum(theo_amp_10_exp, 1e-12))
            theo_db_10_dense = 20 * np.log10(np.maximum(theo_amp_10_dense, 1e-12))
            res_db_10 = exp_db_10 - theo_db_10_exp
            
            # --- Spectral Grid Plot ---
            fig, axs = plt.subplots(2, 2, figsize=(14, 10))
            fig.suptitle(f"Спектральный анализ ({ds})", fontsize=16)
            
            # Top-Left: Linear Angle
            axs[0, 0].plot(angles, exp_amp_05, 'bo', label=f'Exp {f_05:.2f} THz')
            axs[0, 0].plot(angles_dense, theo_amp_05_dense, 'b-', label=f'Theory {f_05:.2f} THz')
            axs[0, 0].plot(angles, exp_amp_10, 'gs', label=f'Exp {f_10:.2f} THz')
            axs[0, 0].plot(angles_dense, theo_amp_10_dense, 'g-', label=f'Theory {f_10:.2f} THz')
            axs[0, 0].set_xlim(-90, 90)
            axs[0, 0].set_title('Угловая зависимость (Лин. шкала)')
            axs[0, 0].set_xlabel('Угол (град)')
            axs[0, 0].set_ylabel('Амплитуда пропускания')
            axs[0, 0].grid(True, alpha=0.5)
            axs[0, 0].legend()
            
            # Bottom-Left: Linear Residuals vs Angle
            axs[1, 0].plot(angles, res_lin_05, 'bo-', label=f'{f_05:.2f} THz')
            axs[1, 0].plot(angles, res_lin_10, 'gs-', label=f'{f_10:.2f} THz')
            axs[1, 0].axhline(0, color='black', linestyle='--', alpha=0.7)
            axs[1, 0].set_xlim(-90, 90)
            axs[1, 0].set_title('Невязка от угла (Лин. шкала)')
            axs[1, 0].set_xlabel('Угол (град)')
            axs[1, 0].set_ylabel('Невязка (Exp - Theory)')
            axs[1, 0].grid(True, alpha=0.5)
            axs[1, 0].legend()
            
            # Top-Right: Log Angle
            axs[0, 1].plot(angles, exp_db_05, 'bo', label=f'Exp {f_05:.2f} THz')
            axs[0, 1].plot(angles_dense, theo_db_05_dense, 'b-', label=f'Theory {f_05:.2f} THz')
            axs[0, 1].plot(angles, exp_db_10, 'gs', label=f'Exp {f_10:.2f} THz')
            axs[0, 1].plot(angles_dense, theo_db_10_dense, 'g-', label=f'Theory {f_10:.2f} THz')
            axs[0, 1].set_xlim(-90, 90)
            axs[0, 1].set_title('Угловая зависимость (Лог. шкала)')
            axs[0, 1].set_xlabel('Угол (град)')
            axs[0, 1].set_ylabel('Пропускание (дБ)')
            axs[0, 1].grid(True, alpha=0.5)
            axs[0, 1].legend()
            
            # Bottom-Right: Log Residuals vs Angle
            axs[1, 1].plot(angles, res_db_05, 'bo-', label=f'{f_05:.2f} THz')
            axs[1, 1].plot(angles, res_db_10, 'gs-', label=f'{f_10:.2f} THz')
            axs[1, 1].axhline(0, color='black', linestyle='--', alpha=0.7)
            axs[1, 1].set_xlim(-90, 90)
            axs[1, 1].set_title('Невязка от угла (Лог. шкала)')
            axs[1, 1].set_xlabel('Угол (град)')
            axs[1, 1].set_ylabel('Невязка в дБ (Exp - Theory)')
            axs[1, 1].grid(True, alpha=0.5)
            axs[1, 1].legend()
            
            plt.tight_layout()
            spec_img_path = images_dir / f"grid_spectral_{ds}.png"
            plt.savefig(spec_img_path)
            plt.close()
            
            # --- 2. Integral Time-Domain Analysis ---
            T_exp_int = []
            bg_spectra_weights = []
            f_basis = None
            
            for a in angles:
                t_s, E_s, t_b, E_b = data_dict[a]
                E_s_clean = E_s - np.mean(E_s)
                E_b_clean = E_b - np.mean(E_b)
                
                E_s_int = np.trapezoid(E_s_clean**2, t_s)
                E_b_int = np.trapezoid(E_b_clean**2, t_b)
                T_exp_int.append(np.sqrt(E_s_int / E_b_int))
                
                f, w = build_spectral_basis(t_b, E_b_clean)
                bg_spectra_weights.append(w)
                if f_basis is None:
                    f_basis = f
                    
            T_exp_int = np.array(T_exp_int)
            weights_avg = np.mean(bg_spectra_weights, axis=0)
            
            alpha_scaled = 0.5 * (LOSS_OPT / 4.343)
            
            # Теоретические значения для экспериментальных углов (для невязок)
            T_theo_int_exp = np.sqrt(np.array([
                theory_T_integral(a, P_OPT, D_OPT, alpha_scaled, GAMMA_OPT, f_basis, weights_avg, OFFSET_OPT)
                for a in angles
            ]))
            
            # Теоретическая кривая для плотной сетки -90..90
            T_theo_int_dense = np.sqrt(np.array([
                theory_T_integral(a, P_OPT, D_OPT, alpha_scaled, GAMMA_OPT, f_basis, weights_avg, OFFSET_OPT)
                for a in angles_dense
            ]))
            
            res_lin_int = T_exp_int - T_theo_int_exp
            
            exp_db_int = 20 * np.log10(np.maximum(T_exp_int, 1e-12))
            theo_db_int_exp = 20 * np.log10(np.maximum(T_theo_int_exp, 1e-12))
            theo_db_int_dense = 20 * np.log10(np.maximum(T_theo_int_dense, 1e-12))
            res_db_int = exp_db_int - theo_db_int_exp
            
            # --- Integral Grid Plot ---
            fig, axs = plt.subplots(2, 2, figsize=(14, 10))
            fig.suptitle(f"Интегральный метод (Time-Domain без DC) - {ds}", fontsize=16)
            
            # Top-Left: Linear Angle
            axs[0, 0].plot(angles, T_exp_int, 'o', color='purple', label='Exp Integral')
            axs[0, 0].plot(angles_dense, T_theo_int_dense, '-', color='purple', label='Theory Integral')
            axs[0, 0].set_xlim(-90, 90)
            axs[0, 0].set_title('Угловая зависимость (Лин. шкала)')
            axs[0, 0].set_xlabel('Угол (град)')
            axs[0, 0].set_ylabel('Амплитуда пропускания')
            axs[0, 0].grid(True, alpha=0.5)
            axs[0, 0].legend()
            
            # Bottom-Left: Linear Residuals vs Angle
            axs[1, 0].plot(angles, res_lin_int, 'o-', color='purple', label='Residual')
            axs[1, 0].axhline(0, color='black', linestyle='--', alpha=0.7)
            axs[1, 0].set_xlim(-90, 90)
            axs[1, 0].set_title('Невязка от угла (Лин. шкала)')
            axs[1, 0].set_xlabel('Угол (град)')
            axs[1, 0].set_ylabel('Невязка (Exp - Theory)')
            axs[1, 0].grid(True, alpha=0.5)
            axs[1, 0].legend()
            
            # Top-Right: Log Angle
            axs[0, 1].plot(angles, exp_db_int, 'o', color='purple', label='Exp Integral')
            axs[0, 1].plot(angles_dense, theo_db_int_dense, '-', color='purple', label='Theory Integral')
            axs[0, 1].set_xlim(-90, 90)
            axs[0, 1].set_title('Угловая зависимость (Лог. шкала)')
            axs[0, 1].set_xlabel('Угол (град)')
            axs[0, 1].set_ylabel('Пропускание (дБ)')
            axs[0, 1].grid(True, alpha=0.5)
            axs[0, 1].legend()
            
            # Bottom-Right: Log Residuals vs Angle
            axs[1, 1].plot(angles, res_db_int, 'o-', color='purple', label='Residual (dB)')
            axs[1, 1].axhline(0, color='black', linestyle='--', alpha=0.7)
            axs[1, 1].set_xlim(-90, 90)
            axs[1, 1].set_title('Невязка от угла (Лог. шкала)')
            axs[1, 1].set_xlabel('Угол (град)')
            axs[1, 1].set_ylabel('Невязка в дБ (Exp - Theory)')
            axs[1, 1].grid(True, alpha=0.5)
            axs[1, 1].legend()
            
            plt.tight_layout()
            int_img_path = images_dir / f"grid_integral_{ds}.png"
            plt.savefig(int_img_path)
            plt.close()
            
            f_report.write(f"## {ds}\n\n")
            f_report.write("### Спектральный анализ (0.5 и 1.0 ТГц)\n")
            f_report.write(f"![Спектральный {ds}](../images/grid_spectral_{ds}.png)\n\n")
            f_report.write("### Интегральный метод (Time-Domain)\n")
            f_report.write(f"![Интегральный {ds}](../images/grid_integral_{ds}.png)\n\n")
            f_report.write("---\n\n")
    print(f"Отчет сгенерирован: {report_path}")

if __name__ == "__main__":
    main()
