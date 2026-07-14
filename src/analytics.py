from pathlib import Path
import matplotlib.pyplot as plt

def generate_report(results: dict, output_path: Path):
    """
    Генерирует Markdown-отчет о результатах оптимизации, сравнивая 1D и 2D методы.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Сводный отчет по оптимизации параметров поляризатора\n\n")
        
        f.write("Данный отчет автоматически сгенерирован конвейером `run_pipeline.py`.\n\n")
        
        f.write("## Сравнение методов 1D (Интегральный) и 2D (Спектрально-угловой)\n\n")
        f.write("| Датасет | Метод | P_eff (мкм) | D_eff (мкм) | Сдвиг (град) | Функция потерь |\n")
        f.write("|---------|-------|-------------|-------------|--------------|----------------|\n")
        
        for ds, ds_results in results.items():
            r1 = ds_results.get('1D')
            if r1:
                f.write(f"| **{ds}** | 1D | {r1['P_eff_um']:.3f} | {r1['D_eff_um']:.3f} | {r1['theta_offset']:.2f} | {r1['fun']:.3e} |\n")
            
            r2 = ds_results.get('2D')
            if r2:
                f.write(f"| | 2D | {r2['P_eff_um']:.3f} | {r2['D_eff_um']:.3f} | {r2['theta_offset']:.2f} | {r2['fun']:.3e} |\n")
                
        f.write("\n## Анализ расхождений (Выбросы и погрешности)\n\n")
        f.write("Сравнение полученных параметров `P_eff` и `D_eff` между методами:\n\n")
        
        for ds, ds_results in results.items():
            if '1D' in ds_results and '2D' in ds_results:
                dp = abs(ds_results['1D']['P_eff_um'] - ds_results['2D']['P_eff_um'])
                dd = abs(ds_results['1D']['D_eff_um'] - ds_results['2D']['D_eff_um'])
                dtheta = abs(ds_results['1D']['theta_offset'] - ds_results['2D']['theta_offset'])
                
                f.write(f"### {ds}\n")
                f.write(f"- ΔP = {dp:.3f} мкм\n")
                f.write(f"- ΔD = {dd:.3f} мкм\n")
                f.write(f"- ΔСдвиг = {dtheta:.2f}°\n")
                if dp > 0.5 or dd > 0.5:
                    f.write("> [!WARNING]\n> Значительное расхождение параметров между методами (более 0.5 мкм). Возможны проблемы с качеством данных или сильным дрейфом лазера в данном датасете.\n\n")
                else:
                    f.write("\n")
                    
        f.write("\n---\n*Отчет сгенерирован автоматически.*\n")
