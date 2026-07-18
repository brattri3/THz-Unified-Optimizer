# -*- coding: utf-8 -*-
"""
scripts/diagnose_40_20.py

Диагностический скрипт для данных test_grid_40_20 (P=40 мкм, D=20 мкм).
Постановка эксперимента: WGP между двумя идеальными плёночными поляризаторами.

Что делает скрипт:
1. Загружает все 25 угловых точек
2. Проверяет симметрию T(+theta) = T(-theta)
3. Строит угловую зависимость пропускания и сравнивает с законом Малюса
4. Рисует мозаику временных трасс
5. Сохраняет диагностический отчёт
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
import scipy.stats as stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from unified_optimizer.data_manager import DataManager
from unified_optimizer import config

# ─── Пути ─────────────────────────────────────────────────────────────────────
DATA_DIR    = Path(config.DATA_DIR)
IMAGES_DIR  = Path(config.BASE_DIR.parent / "docs" / "images")
REPORT_PATH = Path(config.BASE_DIR.parent / "docs" / "artifacts" / "diagnose_40_20_report.md")
IMAGES_DIR.mkdir(parents=True, exist_ok=True)
DATASET     = "test_grid_40_20"

# ─── Загрузка данных ───────────────────────────────────────────────────────────
def load_all_angles(data_dir: Path):
    """Загружает все файлы датасета и возвращает словарь {angle: (time, sig, bg)}."""
    result = {}
    dataset_dir = data_dir / DATASET
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Directory not found: {dataset_dir}")

    for sig_file in sorted(dataset_dir.glob("*_sig.txt")):
        bg_file = Path(str(sig_file).replace("_sig.txt", "_bg.txt"))
        if not bg_file.exists():
            continue
        # Парсим угол из имени файла
        name = sig_file.stem  # e.g. test_grid_40_20_-90deg_rep1_sig
        parts = name.split("_")
        angle_str = None
        for p in parts:
            if "deg" in p:
                angle_str = p.replace("deg", "")
                break
        if angle_str is None:
            continue
        try:
            angle = int(angle_str)
        except ValueError:
            continue

        sig_data = np.loadtxt(sig_file)
        bg_data  = np.loadtxt(bg_file)
        result[angle] = (sig_data[:, 0], sig_data[:, 1], bg_data[:, 1])

    return result


def peak_transmission(sig: np.ndarray, bg: np.ndarray) -> float:
    """Отношение пиковых амплитуд (мера пропускания по времени)."""
    return np.max(np.abs(sig)) / np.max(np.abs(bg))


def malus_law(theta_deg: np.ndarray, t_perp: float, t_par: float) -> np.ndarray:
    """T = (t_perp * cos^2 + t_par * sin^2)^2 — закон Малюса для WGP."""
    th = np.deg2rad(theta_deg)
    return (t_perp * np.cos(th)**2 + t_par * np.sin(th)**2)**2


# ─── Основная логика ──────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("ДИАГНОСТИКА ДАТАСЕТА test_grid_40_20")
    print("Конфигурация: [Film P1] -> [WGP] -> [Film P2]")
    print("=" * 60)

    data = load_all_angles(DATA_DIR)
    angles_sorted = sorted(data.keys())
    print(f"\nЗагружено {len(data)} угловых точек: {angles_sorted}\n")

    # ── 1. Расчёт пиковых пропусканий ─────────────────────────────────────────
    T_all = {a: peak_transmission(data[a][1], data[a][2]) for a in angles_sorted}

    print("Угловые пропускания:")
    print(f"{'Угол':>8}  {'T_peak':>10}")
    print("-" * 22)
    for a in angles_sorted:
        print(f"{a:>7}deg  {T_all[a]:>10.4f}")

    # ── 2. Симметрия T(+theta) vs T(-theta) ───────────────────────────────────
    pos_angles = [a for a in angles_sorted if a > 0]
    print("\n\nПроверка симметрии T(+theta) vs T(-theta):")
    print(f"{'Угол':>8}  {'T(+)':>8}  {'T(-)':>8}  {'Asym%':>8}")
    print("-" * 42)
    asym_list = []
    for a in pos_angles:
        if -a in T_all:
            Tp = T_all[a]
            Tm = T_all[-a]
            mid = (Tp + Tm) / 2
            asym = abs(Tp - Tm) / mid * 100 if mid > 1e-9 else 0.0
            asym_list.append(asym)
            flag = " !!!" if asym > 10 else ""
            print(f"{a:>7}deg  {Tp:>8.4f}  {Tm:>8.4f}  {asym:>7.2f}%{flag}")
    if asym_list:
        print(f"\nСредняя асимметрия: {np.mean(asym_list):.2f}%, "
              f"Максимальная: {np.max(asym_list):.2f}%")
        quality = "ХОРОШО (< 5%)" if np.max(asym_list) < 5 else \
                  "ПРИЕМЛЕМО (< 10%)" if np.max(asym_list) < 10 else \
                  "ВНИМАНИЕ (> 10%)"
        print(f"Оценка симметрии: {quality}")

    # ── 3. Подгонка закона Малюса ──────────────────────────────────────────────
    all_angles_arr = np.array(angles_sorted, dtype=float)
    T_arr = np.array([T_all[a] for a in angles_sorted])

    popt, pcov = curve_fit(malus_law, all_angles_arr, T_arr,
                           p0=[0.99, 0.13], bounds=([0, 0], [2, 1]),
                           maxfev=5000)
    perr = np.sqrt(np.diag(pcov))
    t_perp_fit, t_par_fit = popt
    T_pred = malus_law(all_angles_arr, *popt)
    rmse_malus = np.sqrt(np.mean((T_arr - T_pred)**2))
    r2 = 1 - np.sum((T_arr - T_pred)**2) / np.sum((T_arr - np.mean(T_arr))**2)

    print(f"\n\nПодгонка закона Малюса: T = (t_perp*cos^2 + t_par*sin^2)^2")
    print(f"  t_perp = {t_perp_fit:.5f} +/- {perr[0]:.5f}")
    print(f"  t_par  = {t_par_fit:.5f} +/- {perr[1]:.5f}")
    print(f"  t_perp^2 = {t_perp_fit**2:.5f}  (интенсивность на 0 deg)")
    print(f"  t_par^2  = {t_par_fit**2:.5f}  (интенсивность на 90 deg)")
    print(f"  Extinction ratio = {t_perp_fit**2 / t_par_fit**2:.1f}:1")
    print(f"  RMSE = {rmse_malus:.5f},  R^2 = {r2:.6f}")

    # ── 4. Графики ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Diagnostics: {DATASET} (P=40 um, D=20 um)\n"
                 "[Film P1] -> [WGP] -> [Film P2]", fontsize=13, fontweight='bold')

    # --- 4a. Угловая зависимость T(theta) ---
    ax = axes[0, 0]
    ax.scatter(all_angles_arr, T_arr, color='steelblue', zorder=5,
               label='Experimental data', s=60)
    theta_fine = np.linspace(-105, 105, 400)
    ax.plot(theta_fine, malus_law(theta_fine, *popt), 'r-', linewidth=2,
            label=f'Malus fit: t_perp={t_perp_fit:.3f}, t_par={t_par_fit:.3f}')
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("Peak Transmission T")
    ax.set_title("Angular Dependence T(theta)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

    # --- 4b. Симметрия T(+θ) vs T(-θ) ---
    ax = axes[0, 1]
    sym_pos = [T_all[a] for a in pos_angles if -a in T_all]
    sym_neg = [T_all[-a] for a in pos_angles if -a in T_all]
    ax.scatter(sym_pos, sym_neg, color='darkorange', s=60)
    lim = max(max(sym_pos), max(sym_neg)) * 1.05
    ax.plot([0, lim], [0, lim], 'k--', linewidth=1, label='Perfect symmetry')
    ax.set_xlabel("T(+theta)")
    ax.set_ylabel("T(-theta)")
    ax.set_title("Symmetry Check: T(+theta) vs T(-theta)")
    ax.legend()
    ax.grid(True, alpha=0.4)

    # --- 4c. Мозаика временных трасс (только положительные углы) ---
    ax = axes[1, 0]
    pos_only = sorted([a for a in angles_sorted if a >= 0])
    cmap = plt.cm.viridis(np.linspace(0, 1, len(pos_only)))
    for i, a in enumerate(pos_only):
        t, sig, bg = data[a]
        T_trace = sig / (np.max(np.abs(bg)) + 1e-15)
        ax.plot(t, T_trace, color=cmap[i], linewidth=0.7, alpha=0.85,
                label=f'{a}°' if a % 30 == 0 else None)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Normalized Amplitude")
    ax.set_title("Time-domain traces (positive angles)")
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([t[0], t[-1]])

    # --- 4d. Невязки закона Малюса ---
    ax = axes[1, 1]
    residuals = T_arr - T_pred
    ax.bar(all_angles_arr, residuals, color='slateblue', alpha=0.7, width=4)
    ax.axhline(0, color='k', linewidth=1)
    ax.axhline( 2*rmse_malus, color='r', linestyle='--', linewidth=1, label=f'+/-2*RMSE ({2*rmse_malus:.4f})')
    ax.axhline(-2*rmse_malus, color='r', linestyle='--', linewidth=1)
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("Residual T_exp - T_model")
    ax.set_title("Residuals of Malus Law Fit")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    fig_path = IMAGES_DIR / "diagnose_40_20.png"
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f"\nГрафик сохранен: {fig_path}")

    # ── 5. Сохранение отчёта ─────────────────────────────────────────────────
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write(f"# Диагностический отчёт: {DATASET}\n\n")
        f.write(f"**Дата**: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write("**Конфигурация эксперимента**: `[Ideal Film P1] → [WGP (P=40, D=20 мкм)] → [Ideal Film P2]`\n\n")
        f.write(f"**Количество угловых точек**: {len(data)}\n\n")
        f.write("## Подгонка закона Малюса\n\n")
        f.write(f"| Параметр | Значение |\n|---|---|\n")
        f.write(f"| t_perp | {t_perp_fit:.5f} ± {perr[0]:.5f} |\n")
        f.write(f"| t_par | {t_par_fit:.5f} ± {perr[1]:.5f} |\n")
        f.write(f"| Extinction ratio | {t_perp_fit**2/t_par_fit**2:.1f}:1 |\n")
        f.write(f"| RMSE | {rmse_malus:.5f} |\n")
        f.write(f"| R² | {r2:.6f} |\n\n")
        f.write(f"## Симметрия\n\n")
        if asym_list:
            f.write(f"Средняя асимметрия T(+θ)/T(−θ): **{np.mean(asym_list):.2f}%**  \n")
            f.write(f"Максимальная асимметрия: **{np.max(asym_list):.2f}%**  \n")
            f.write(f"Оценка: **{quality}**\n\n")
        f.write(f"![Diagnostics plot]({fig_path})\n")

    print(f"Отчёт сохранен: {REPORT_PATH}")
    print("\nДиагностика завершена.")


if __name__ == "__main__":
    main()
