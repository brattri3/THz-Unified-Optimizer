import os
import argparse
from pathlib import Path

def generate_templates(dataset_name: str, out_dir: str = "data_pool"):
    """
    Генерирует пустые текстовые файлы с правильными именами для новой серии ТГц измерений.
    """
    # Базовые углы: от -90 до 90 с шагом 10
    base_angles = list(range(-90, 100, 10))
    
    # Дополнительные углы для скрещенного состояния (около -90 и +90 с шагом 2 градуса)
    fine_angles_minus = list(range(-100, -78, 2))
    fine_angles_plus = list(range(80, 102, 2))
    
    # Объединяем, удаляем дубликаты и сортируем
    all_angles = sorted(list(set(base_angles + fine_angles_minus + fine_angles_plus)))
    
    # Ограничиваем физическим диапазоном от -180 до 180
    all_angles = [a for a in all_angles if -180 <= a <= 180]
    
    # Создаем папку под серию
    target_dir = Path(out_dir) / dataset_name
    target_dir.mkdir(parents=True, exist_ok=True)
    
    created_count = 0
    for angle in all_angles:
        # Для каждого угла создаем файл сигнала (sig) и фона (bg)
        sig_file = target_dir / f"{dataset_name}_{angle}deg_rep1_sig.txt"
        bg_file = target_dir / f"{dataset_name}_{angle}deg_rep1_bg.txt"
        
        # Создаем пустые файлы
        sig_file.touch()
        bg_file.touch()
        created_count += 2
        
    print(f"Успешно создано {created_count} пустых файлов-шаблонов в директории:\n{target_dir.absolute()}")
    print("Вам осталось только скопировать в них столбцы (Время, Сигнал) с THz-TDS установки.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Генератор пустых файлов-шаблонов для ТГц измерений.")
    parser.add_argument("dataset_name", type=str, nargs="?", default="grid_40_20",
                        help="Имя новой серии (по умолчанию: grid_40_20)")
    parser.add_argument("--out", type=str, default="data_pool",
                        help="Папка для сохранения (по умолчанию: data_pool)")
    
    args = parser.parse_args()
    generate_templates(args.dataset_name, args.out)
