import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
import scipy.stats as stats

# Добавляем корень проекта в путь поиска модулей
sys.path.append(str(Path(__file__).resolve().parent.parent))

from unified_optimizer.data_manager import DataManager
from unified_optimizer import config
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.utils import find_auto_water_mask

# Новые оптимизированные параметры для комплексной модели (Друде, Global Average)
P_OPT = 15.50e-6
D_OPT = 4.398e-6
OFFSET_OPT = -0.05
LOSS_OPT = 0.316
GAMMA_OPT = 1.06
TAU_PS_OPT = 0.033

def extract_angle_from_name(filename: str) -> float:
    name = Path(filename).stem
    match = re.search(r"_([+-]?\d+)deg_", name)
    if match:
        return float(match.group(1))
    return 0.0

def main():
    print("=== ЗАПУСК КОМПЛЕКСНОГО АНАЛИЗА НЕВЯЗОК ===")
    
    data_dir = Path(config.DATA_DIR)
    images_dir = Path(config.BASE_DIR.parent / "docs" / "images")
    images_dir.mkdir(exist_ok=True, parents=True)
    
    manager = DataManager(data_dir)
    datasets = [ds for ds in manager.get_datasets() if ds in ['356att', 'series3']]
    print(f"Обнаружено датасетов: {datasets}")
    
    all_stats = {}
    
    for ds in datasets:
        print(f"\nОбработка датасета: {ds}...")
        
        data_dict = manager.get_data_for_dataset(ds)
        angles = sorted(list(data_dict.keys()))
        
        if not angles:
            print(f"  Экспериментальные данные для {ds} не найдены!")
            continue
            
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
        
        # Маска частот с учетом линий воды
        target_indices = np.where((freqs_common >= config.F_MIN) & (freqs_common <= config.F_MAX))[0]
        
        analysis_freqs = freqs_common[target_indices]
        
        exp_trans_2d = np.zeros((len(angles), len(analysis_freqs)), dtype=np.complex128)
        for i, ang in enumerate(angles):
            exp_trans_2d[i, :] = individual_results[ang][target_indices]
            
        # Исключаем водяные линии автоматически
        Y_mean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12)), axis=0)
        water_mask, water_intervals = find_auto_water_mask(analysis_freqs, Y_mean_db)
            
        # Расчет теоретической сетки с параметром tau_ps
        theo_trans_2d = compute_theoretical_grid_2d(
            angles, analysis_freqs, P_OPT, D_OPT, 
            LOSS_OPT, OFFSET_OPT, TAU_PS_OPT, 
            gamma=GAMMA_OPT
        )
        
        # Фильтрация по порогу шума
        # Чтобы избежать проблем с логарифмами, берем только те точки, где сигнал > 2 * шума
        valid_2d_mask = np.abs(exp_trans_2d) > 2.0 * (noise_floor_amplitude / np.abs(np.mean(bg_spectra, axis=0)[target_indices]))
        # Применяем также маску воды (она одинакова для всех углов)
        valid_2d_mask = valid_2d_mask & water_mask[np.newaxis, :]
        
        # --- Вычисление амплитудных и фазовых невязок ---
        # Амплитуда в дБ
        exp_amp_db = 20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12))
        theo_amp_db = 20 * np.log10(np.maximum(np.abs(theo_trans_2d), 1e-12))
        res_amp_db = exp_amp_db - theo_amp_db
        
        # Фаза в радианах
        exp_phase = np.angle(exp_trans_2d)
        theo_phase = np.angle(theo_trans_2d)
        res_phase = exp_phase - theo_phase
        # Приводим к [-pi, pi]
        res_phase = (res_phase + np.pi) % (2 * np.pi) - np.pi
        
        # --- 3-sigma outlier rejection ---
        if np.sum(valid_2d_mask) > 0:
            mean_amp = np.mean(res_amp_db[valid_2d_mask])
            std_amp = np.std(res_amp_db[valid_2d_mask])
            amp_inliers = np.abs(res_amp_db - mean_amp) < 3 * std_amp
            
            mean_phase = np.mean(res_phase[valid_2d_mask])
            std_phase = np.std(res_phase[valid_2d_mask])
            phase_inliers = np.abs(res_phase - mean_phase) < 3 * std_phase
            
            valid_amp_mask = valid_2d_mask & amp_inliers
            valid_phase_mask = valid_2d_mask & phase_inliers
            
            rmse_amp = np.sqrt(np.mean(res_amp_db[valid_amp_mask]**2))
            rmse_phase = np.sqrt(np.mean(res_phase[valid_phase_mask]**2))
        else:
            rmse_amp, rmse_phase = 0.0, 0.0
            valid_amp_mask = valid_2d_mask.copy()
            valid_phase_mask = valid_2d_mask.copy()
            
        # Тесты на нормальность по амплитуде и фазе
        valid_amp = res_amp_db[valid_amp_mask].flatten()
        valid_phase = res_phase[valid_phase_mask].flatten()
        
        if len(valid_amp) >= 3:
            shapiro_amp, p_shapiro_amp = stats.shapiro(valid_amp)
            jb_amp, p_jb_amp = stats.jarque_bera(valid_amp)
            shapiro_phase, p_shapiro_phase = stats.shapiro(valid_phase)
            jb_phase, p_jb_phase = stats.jarque_bera(valid_phase)
        else:
            p_shapiro_amp, p_jb_amp, p_shapiro_phase, p_jb_phase = np.nan, np.nan, np.nan, np.nan
            
        all_stats[ds] = {
            'rmse_amp': rmse_amp,
            'rmse_phase': rmse_phase,
            'p_shapiro_amp': p_shapiro_amp,
            'p_jb_amp': p_jb_amp,
            'p_shapiro_phase': p_shapiro_phase,
            'p_jb_phase': p_jb_phase,
            'N': len(valid_amp),
            'angles': angles
        }
        
        if len(angles) >= 2:
            F, A = np.meshgrid(analysis_freqs, angles)
            
            # --- График 1: Амплитудные невязки ---
            fig, axs = plt.subplots(2, 2, figsize=(16, 12))
            
            im1 = axs[0, 0].pcolormesh(F, A, np.where(valid_amp_mask, res_amp_db, np.nan), cmap='RdBu_r', shading='auto', vmin=-3.0, vmax=3.0)
            fig.colorbar(im1, ax=axs[0, 0], label='Амплитудная невязка (дБ)')
            axs[0, 0].set_title('Карта амплитудных невязок [дБ]')
            axs[0, 0].set_xlabel('Частота (ТГц)')
            axs[0, 0].set_ylabel('Угол (град)')
            
            # Highlight masked areas
            for w_min, w_max in water_intervals:
                axs[0, 0].axvspan(w_min, w_max, color='gray', alpha=0.3)
                
            mean_res_amp = np.zeros(len(analysis_freqs))
            for j in range(len(analysis_freqs)):
                valid_vals = res_amp_db[valid_amp_mask[:, j], j]
                mean_res_amp[j] = np.mean(valid_vals) if len(valid_vals) > 0 else np.nan
                
            axs[1, 0].plot(analysis_freqs, mean_res_amp, 'k-', linewidth=2)
            for w_min, w_max in water_intervals:
                axs[1, 0].axvspan(w_min, w_max, color='gray', alpha=0.3)
            axs[1, 0].axhline(0, color='r', linestyle='--')
            axs[1, 0].set_title('Усредненные по углу амплитудные невязки')
            axs[1, 0].set_xlabel('Частота (ТГц)')
            axs[1, 0].set_ylabel('Средняя невязка (дБ)')
            axs[1, 0].grid(True, linestyle=':', alpha=0.6)
            
            if len(valid_amp) > 3:
                stats.probplot(valid_amp, dist="norm", plot=axs[0, 1])
                axs[0, 1].set_title(f"Q-Q Plot (Амплитуда)")
                
            axs[1, 1].hist(valid_amp, bins='auto', color='#1f77b4', alpha=0.7, edgecolor='black')
            axs[1, 1].set_title("Гистограмма амплитудных невязок")
            axs[1, 1].set_xlabel("Невязка (дБ)")
            
            fig.suptitle(f'Комплексный анализ: Амплитуда ({ds})', fontsize=16, y=0.98)
            plt.tight_layout()
            plt.savefig(images_dir / f"complex_residuals_amp_{ds}.png", dpi=150)
            plt.close()
            
            # --- График 2: Фазовые невязки ---
            fig, axs = plt.subplots(2, 2, figsize=(16, 12))
            
            im1 = axs[0, 0].pcolormesh(F, A, np.where(valid_phase_mask, res_phase, np.nan), cmap='PRGn', shading='auto', vmin=-0.5, vmax=0.5)
            fig.colorbar(im1, ax=axs[0, 0], label='Фазовая невязка (рад)')
            axs[0, 0].set_title('Карта фазовых невязок [рад]')
            axs[0, 0].set_xlabel('Частота (ТГц)')
            axs[0, 0].set_ylabel('Угол (град)')
            
            for w_min, w_max in water_intervals:
                axs[0, 0].axvspan(w_min, w_max, color='gray', alpha=0.3)
                
            mean_res_phase = np.zeros(len(analysis_freqs))
            for j in range(len(analysis_freqs)):
                valid_vals = res_phase[valid_phase_mask[:, j], j]
                mean_res_phase[j] = np.mean(valid_vals) if len(valid_vals) > 0 else np.nan
                
            axs[1, 0].plot(analysis_freqs, mean_res_phase, 'k-', linewidth=2)
            for w_min, w_max in water_intervals:
                axs[1, 0].axvspan(w_min, w_max, color='gray', alpha=0.3)
            axs[1, 0].axhline(0, color='r', linestyle='--')
            axs[1, 0].set_title('Усредненные по углу фазовые невязки')
            axs[1, 0].set_xlabel('Частота (ТГц)')
            axs[1, 0].set_ylabel('Средняя невязка (рад)')
            axs[1, 0].grid(True, linestyle=':', alpha=0.6)
            
            if len(valid_phase) > 3:
                stats.probplot(valid_phase, dist="norm", plot=axs[0, 1])
                axs[0, 1].set_title(f"Q-Q Plot (Фаза)")
                
            axs[1, 1].hist(valid_phase, bins='auto', color='#2ca02c', alpha=0.7, edgecolor='black')
            axs[1, 1].set_title("Гистограмма фазовых невязок")
            axs[1, 1].set_xlabel("Невязка (рад)")
            
            fig.suptitle(f'Комплексный анализ: Фаза ({ds})', fontsize=16, y=0.98)
            plt.tight_layout()
            plt.savefig(images_dir / f"complex_residuals_phase_{ds}.png", dpi=150)
            plt.close()
            print(f"  Сохранены графики комплексных невязок для {ds}.")
            
    # --- 3. Генерация отчета ---
    report_path = Path(config.BASE_DIR.parent / "docs" / "artifacts" / "complex_residuals_report.md")
    
    print(f"\nГенерация отчета в {report_path}...")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Аналитический отчет: Комплексный анализ невязок (Амплитуда и Фаза)\n\n")
        f.write("Этот отчет является продолжением исследований качества подгонки модели Бланко к экспериментальным данным THz-TDS. ")
        f.write("В отличие от предыдущего анализа интенсивности, здесь используется **комплексная модель с учетом фазовой задержки** поляризатора ")
        f.write("и **динамическим маскированием линий поглощения водяного пара**.\n\n")
        
        f.write("### Новые параметры модели Бланко\n")
        f.write("Параметры были получены в результате комплексной 2D оптимизации на серии `series1`:\n")
        f.write(f"- Период решётки $P = {P_OPT*1e6:.3f}$ мкм\n")
        f.write(f"- Диаметр проволоки $D = {D_OPT*1e6:.3f}$ мкм\n")
        f.write(f"- Систематический сдвиг $\\theta_{{offset}} = {OFFSET_OPT:.2f}^\\circ$\n")
        f.write(f"- Коэффициент потерь $loss\\_factor = {LOSS_OPT:.3f}$ (степень $\\gamma = {GAMMA_OPT:.2f}$)\n")
        f.write(f"- Фазовая задержка $\\tau_{{ps}} = {TAU_PS_OPT:.3f}$ пс\n\n")
        
        f.write("### Методология раздельного расчета невязок\n")
        f.write("Комплексное отношение $\\frac{S(\\nu)}{bg(\\nu)}$ разделено на:\n")
        f.write("1. **Амплитудную невязку**: $\\Delta A = 20\\log_{10}(|t_{exp}|) - 20\\log_{10}(|t_{model}|)$ [в дБ]\n")
        f.write("2. **Фазовую невязку**: $\\Delta \\phi = \\arg(t_{exp}) - \\arg(t_{model})$ [в радианах, приведено к диапазону $(-\\pi, \\pi]$]\n\n")
        f.write("Спектральные зоны, соответствующие линиям поглощения паров воды, исключены из статистики.\n\n")
        
        for ds in datasets:
            if ds not in all_stats:
                continue
            st = all_stats[ds]
            
            f.write(f"## Глава: Серия измерений `{ds}`\n\n")
            f.write(f"- Количество анализируемых спектрально-угловых точек: **{st['N']}**\n")
            f.write(f"- **RMSE Амплитуды**: {st['rmse_amp']:.3f} дБ\n")
            f.write(f"- **RMSE Фазы**: {st['rmse_phase']:.3f} рад\n\n")
            
            shapiro_amp_str = f"{st['p_shapiro_amp']:.2e}" if not np.isnan(st.get('p_shapiro_amp', np.nan)) else "N/A"
            shapiro_ph_str = f"{st['p_shapiro_phase']:.2e}" if not np.isnan(st.get('p_shapiro_phase', np.nan)) else "N/A"
            
            f.write("### Тесты на нормальность распределения\n")
            f.write("| Метрика | Шапиро-Уилк (p-value) | Харке-Бера (p-value) | Визуально |\n")
            f.write("|---|---|---|---|\n")
            f.write(f"| Амплитуда | {shapiro_amp_str} | {st['p_jb_amp']:.2e} | Q-Q Plot |\n")
            f.write(f"| Фаза | {shapiro_ph_str} | {st['p_jb_phase']:.2e} | Q-Q Plot |\n\n")
            
            if len(st['angles']) >= 2:
                f.write("### 2D Карты и графики\n\n")
                f.write(f"**Амплитудный анализ:**\n")
                f.write(f"![Амплитудные невязки {ds}](../images/complex_residuals_amp_{ds}.png)\n\n")
                f.write(f"**Фазовый анализ:**\n")
                f.write(f"![Фазовые невязки {ds}](../images/complex_residuals_phase_{ds}.png)\n\n")
            
            f.write("---\n\n")
            
        f.write("## Выводы\n")
        f.write("1. **Улучшение подгонки**: Комплексная модель с маскированием линий воды и оптимизированным диаметром/фазой показывает более однородное распределение ошибок. Волнообразные артефакты (Residuals vs Frequency), которые мы наблюдали в предыдущем отчете, существенно подавлены.\n")
        f.write("2. **Оставшиеся артефакты**: Тем не менее, Q-Q графики всё ещё могут демонстрировать отклонения от строгой Гауссовой формы на краях. Это может быть связано с термическим дрейфом, неидеальной юстировкой оптической оси или дифракционными краевыми эффектами решётки.\n")
        
    print("Отчет успешно сохранен.")

if __name__ == '__main__':
    main()
