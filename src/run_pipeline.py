from pathlib import Path
from unified_optimizer import config
from unified_optimizer.data_manager import DataManager
from unified_optimizer.optimizer_1d import optimize_1d_integral
from unified_optimizer.optimizer_2d import optimize_2d_spectral
from unified_optimizer.analytics import generate_report

def main():
    print("Инициализация DataManager...")
    manager = DataManager(config.DATA_DIR)
    datasets = manager.get_datasets()
    
    if not datasets:
        print(f"Ошибка: Не найдено датасетов в {config.DATA_DIR}")
        return
        
    print(f"Найдено {len(datasets)} датасетов: {datasets}")
    
    results = {}
    
    # Обработка индивидуальных датасетов
    for ds in datasets:
        print(f"\n--- Обработка датасета: {ds} ---")
        data_dict = manager.get_data_for_dataset(ds)
        angles = sorted(list(data_dict.keys()))
        print(f"Загружено {len(angles)} углов: {angles}")
        
        print("  Запуск 1D Интегральной оптимизации...")
        res_1d = optimize_1d_integral(data_dict)
        print(f"    P = {res_1d['P_eff_um']:.3f} мкм, D = {res_1d['D_eff_um']:.3f} мкм, "
              f"Сдвиг = {res_1d['theta_offset']:.2f}°")
              
        print("  Запуск 2D Спектральной оптимизации...")
        res_2d = optimize_2d_spectral(data_dict)
        print(f"    P = {res_2d['P_eff_um']:.3f} мкм, D = {res_2d['D_eff_um']:.3f} мкм, "
              f"Сдвиг = {res_2d['theta_offset']:.2f}°")
              
        results[ds] = {
            '1D': res_1d,
            '2D': res_2d
        }
        
    # Обработка глобального усреднения
    print("\n--- Обработка глобального усреднения (Global Average) ---")
    global_dict = manager.get_global_average()
    print("  Запуск 1D Интегральной оптимизации...")
    res_1d_glob = optimize_1d_integral(global_dict)
    
    print("  Запуск 2D Спектральной оптимизации...")
    res_2d_glob = optimize_2d_spectral(global_dict)
    
    results['Global_Average'] = {
        '1D': res_1d_glob,
        '2D': res_2d_glob
    }
    
    # Генерация отчета
    report_path = config.RESULTS_DIR / "optimization_report.md"
    generate_report(results, report_path)
    print(f"\n✅ Отчет успешно сгенерирован: {report_path}")

if __name__ == "__main__":
    main()
