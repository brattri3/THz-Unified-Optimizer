import os
from pathlib import Path
import numpy as np
import re
from typing import Dict, List, Tuple

def parse_filename(filename: str) -> tuple:
    """
    Парсит имя нормализованного файла вида {dataset}_{angle}deg_rep{N}_{type}.txt
    Например: series1_-10deg_rep1_sig.txt
    Возвращает (dataset_name, angle_deg, rep, type_)
    """
    name = Path(filename).stem
    match = re.match(r"^([^_]+)_([+-]?\d+)deg_rep(\d+)_(sig|bg)$", name)
    if not match:
        raise ValueError(f"Неверный формат имени файла: {filename}")
    dataset_name = match.group(1)
    angle_deg = float(match.group(2))
    rep = int(match.group(3))
    type_ = match.group(4)
    return dataset_name, angle_deg, rep, type_

def load_tds(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Загружает файл ТГц-данных, вычитает постоянную составляющую."""
    raw = np.loadtxt(filepath)
    t, E = raw[:, 0], raw[:, 1]
    E -= np.mean(E[:50])
    return t, E

class DataManager:
    def __init__(self, data_dir: Path):
        self.data_dir = Path(data_dir)
        self.raw_data = {}  # (dataset, angle, rep) -> {'sig': (t, E), 'bg': (t, E)}
        self.scan_directory()

    def scan_directory(self):
        for f in self.data_dir.glob("*.txt"):
            try:
                dataset_name, angle, rep, type_ = parse_filename(f.name)
            except ValueError:
                continue
            
            key = (dataset_name, angle, rep)
            if key not in self.raw_data:
                self.raw_data[key] = {}
            
            t, E = load_tds(f)
            self.raw_data[key][type_] = (t, E)

    def get_datasets(self) -> List[str]:
        return sorted(list(set(k[0] for k in self.raw_data.keys())))

    def get_data_for_dataset(self, dataset: str) -> Dict[float, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """
        Возвращает данные для конкретного датасета (серии).
        Возврат: dict[angle] = (t_s, E_s, t_b, E_b)
        Если для одного угла есть несколько повторений в одном датасете, они усредняются.
        """
        grouped = {}
        for (ds, angle, rep), data in self.raw_data.items():
            if ds == dataset and 'sig' in data and 'bg' in data:
                if angle not in grouped:
                    grouped[angle] = []
                grouped[angle].append(data)
                
        result = {}
        for angle, entries in grouped.items():
            if len(entries) == 1:
                t_s, E_s = entries[0]['sig']
                t_b, E_b = entries[0]['bg']
            else:
                # Усредняем спектры
                t_s = entries[0]['sig'][0]
                t_b = entries[0]['bg'][0]
                E_s = np.mean([e['sig'][1] for e in entries], axis=0)
                E_b = np.mean([e['bg'][1] for e in entries], axis=0)
            result[angle] = (t_s, E_s, t_b, E_b)
            
        return result

    def get_global_average(self) -> Dict[float, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """
        Усредняет сырые спектры всех доступных серий для каждого уникального угла.
        """
        grouped = {}
        for (ds, angle, rep), data in self.raw_data.items():
            if 'sig' in data and 'bg' in data:
                if angle not in grouped:
                    grouped[angle] = []
                grouped[angle].append(data)
                
        result = {}
        for angle, entries in grouped.items():
            t_s = entries[0]['sig'][0]
            t_b = entries[0]['bg'][0]
            E_s = np.mean([e['sig'][1] for e in entries], axis=0)
            E_b = np.mean([e['bg'][1] for e in entries], axis=0)
            result[angle] = (t_s, E_s, t_b, E_b)
            
        return result
