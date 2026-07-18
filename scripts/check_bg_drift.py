# -*- coding: utf-8 -*-
"""
scripts/check_bg_drift.py

Проверка стабильности / дрейфа лазера по фоновым файлам (*_bg.txt) датасета test_grid_40_20.

Алгоритм:
1. Для каждого *_bg.txt читается время создания/модификации файла (mtime).
2. Вычитается DC-составляющая (среднее по первому участку до основного импульса).
3. Вычисляется интеграл квадрата напряженности поля: Energy = int E^2(t) dt.
4. Строится график зависмости энергии фонового сигнала от времени съемки.
5. Рассчитывается относительное стандартное отклонение (RSD) и дрейф в %.
"""
import os
import sys
from pathlib import Path
import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from unified_optimizer import config

DATASET_DIR = Path(config.DATA_DIR) / "test_grid_40_20"
IMAGES_DIR = Path(config.BASE_DIR.parent / "docs" / "images")
ARTIFACTS_DIR = Path(config.BASE_DIR.parent / "docs" / "artifacts")

IMAGES_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

def main():
    bg_files = sorted(DATASET_DIR.glob("*_bg.txt"))
    if not bg_files:
        print("Файлы bg не найдены!")
        return

    records = []

    for f in bg_files:
        mtime = f.stat().st_mtime
        dt_obj = datetime.datetime.fromtimestamp(mtime)

        # Выделяем угол из имени для подписи
        name = f.stem
        angle_str = None
        for p in name.split("_"):
            if "deg" in p:
                angle_str = p.replace("deg", "")
                break
        angle = int(angle_str) if angle_str else 0

        raw = np.loadtxt(f)
        t = raw[:, 0]  # ps
        E = raw[:, 1]

        # 1. Вычитание DC составляющей (по первым 50 точкам)
        dc_offset = np.mean(E[:50])
        E_clean = E - dc_offset

        # 2. Интегрирование квадрата амплитуды по времени (методом трапеций)
        # int E^2 dt в условных единицах
        energy = np.trapezoid(E_clean**2, t)

        records.append({
            'file': f.name,
            'mtime': mtime,
            'dt': dt_obj,
            'angle': angle,
            'dc_offset': dc_offset,
            'peak_amp': np.max(np.abs(E_clean)),
            'energy': energy
        })

    # Сортируем записи по реальному времени съёмки (mtime)
    records.sort(key=lambda x: x['mtime'])

    start_time = records[0]['mtime']
    times_min = np.array([(r['mtime'] - start_time) / 60.0 for r in records])
    energies = np.array([r['energy'] for r in records])
    angles = [r['angle'] for r in records]
    peaks = np.array([r['peak_amp'] for r in records])

    # Статистика дрейфа энергии
    e_mean = np.mean(energies)
    e_std = np.std(energies)
    rsd_energy = (e_std / e_mean) * 100.0
    total_drift_energy = (energies[-1] - energies[0]) / energies[0] * 100.0
    max_dev_energy = (np.max(energies) - np.min(energies)) / e_mean * 100.0

    # Статистика дрейфа пиковой амплитуды
    p_mean = np.mean(peaks)
    rsd_peak = (np.std(peaks) / p_mean) * 100.0

    print("=" * 65)
    print("АНАЛИЗ ДРЕЙФА ИНТЕГРАЛЬНОЙ ЭНЕРГИИ И АМПЛИТУДЫ ФОНА (BG)")
    print("=" * 65)
    print(f"Количество измерявшихся bg-файлов: {len(records)}")
    print(f"Общая длительность эксперимента:  {times_min[-1]:.1f} мин")
    print(f"Время начала съёмки:             {records[0]['dt'].strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Время окончания съёмки:          {records[-1]['dt'].strftime('%Y-%m-%d %H:%M:%S')}\n")

    print(f"Средняя интегральная энергия int E^2 dt: {e_mean:.5f}")
    print(f"Относительное стд. отклонение (RSD E):   {rsd_energy:.2f}%")
    print(f"Полный дрейф энергии (start -> end):     {total_drift_energy:+.2f}%")
    print(f"Максимальный разброс энергии (max-min):  {max_dev_energy:.2f}%\n")

    print(f"Средняя пиковая амплитуда max|E|:        {p_mean:.5f}")
    print(f"Относительное стд. отклонение (RSD Peak): {rsd_peak:.2f}%\n")

    print(f"{'N':>2} | {'Time':>8} | {'Angle':>6} | {'Peak |E|':>10} | {'Integr E^2':>12} | {'dE from avg':>14}")
    print("-" * 65)
    for i, r in enumerate(records):
        t_m = (r['mtime'] - start_time) / 60.0
        dev = (r['energy'] - e_mean) / e_mean * 100.0
        print(f"{i+1:2d} | {t_m:6.1f}m | {r['angle']:5d}° | {r['peak_amp']:10.4f} | {r['energy']:12.5f} | {dev:+13.2f}%")

    # ── Графики ───────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Laser Stability & Drift Analysis (Background files)", fontsize=13, fontweight='bold')

    # 1. Интеграл квадрата поля int E^2 dt
    ax1.plot(times_min, energies, 'o-', color='crimson', linewidth=1.5, markersize=6, label='Integral E^2 dt')
    ax1.axhline(e_mean, color='k', linestyle='--', label=f'Mean = {e_mean:.4f}')
    ax1.fill_between(times_min, e_mean - e_std, e_mean + e_std, color='crimson', alpha=0.15, label=f'+/-1 Std (RSD={rsd_energy:.2f}%)')
    ax1.set_ylabel("Pulse Energy (arb. u.)")
    ax1.set_title("Time-Integrated Pulse Energy vs Time")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.4)

    # 2. Пиковая амплитуда max|E|
    ax2.plot(times_min, peaks, 's-', color='navy', linewidth=1.5, markersize=5, label='Peak Amplitude |E|')
    ax2.axhline(p_mean, color='k', linestyle='--', label=f'Mean Peak = {p_mean:.4f}')
    ax2.set_xlabel("Elapsed Time (minutes)")
    ax2.set_ylabel("Peak Amplitude (arb. u.)")
    ax2.set_title("Peak Background Amplitude vs Time")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.4)

    # Подписи углов над точками
    for t_m, e_v, a in zip(times_min, energies, angles):
        ax1.annotate(f"{a}°", (t_m, e_v), textcoords="offset points", xytext=(0,6), ha='center', fontsize=7, alpha=0.8)

    plt.tight_layout()
    fig_path = IMAGES_DIR / "laser_drift_bg_40_20.png"
    plt.savefig(fig_path, dpi=150)
    plt.close()

    # Сохраняем текстовый артефакт
    report_path = ARTIFACTS_DIR / "laser_drift_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Анализ стабильности и дрейфа лазера ТГц-спектрометра\n\n")
        f.write(f"**Дата анализа**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**Датасет**: `test_grid_40_20` (всего {len(records)} BG-файлов)\n")
        f.write(f"**Длительность эксперимента**: {times_min[-1]:.1f} минут\n\n")
        f.write("## Статистика стабильности фона\n\n")
        f.write(f"| Метрика | Значение |\n|---|---|\n")
        f.write(f"| Средняя энергия $\\int E^2 dt$ | {e_mean:.5f} |\n")
        f.write(f"| **Относительный разброс (RSD энергии)** | **{rsd_energy:.2f}%** |\n")
        f.write(f"| Суммарный дрейф (начало → конец) | {total_drift_energy:+.2f}% |\n")
        f.write(f"| Максимальное отклонение (max−min) | {max_dev_energy:.2f}% |\n")
        f.write(f"| **RSD пиковой амплитуды** | **{rsd_peak:.2f}%** |\n\n")
        f.write(f"![Laser Drift Plot]({fig_path})\n")

    print(f"\nГрафик сохранен: {fig_path}")
    print(f"Отчет сохранен: {report_path}")

if __name__ == "__main__":
    main()
