# -*- coding: utf-8 -*-
"""Графики для иллюстрированного отчёта A19 (микрофото WGP).

Строит рисунки ТОЛЬКО из готовых артефактов `a18_*.json` — никаких новых расчётов
по изображениям здесь нет. Автор чисел — зона A (заходы A4/A5); оформление — ORCH
по прямой просьбе владельца 2026-08-10.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe research/results/a18_microscopy/figs_report/make_report_figs.py

Палитра — валидированная категориальная (dataviz): слоты blue/orange/aqua/yellow,
проверены `validate_palette.js` (light, adjacent и all-pairs для первых трёх).
Контраст aqua/yellow к фону ниже 3:1 -> везде есть прямые подписи или легенда.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

REPO = Path(__file__).resolve().parents[4]
RES = REPO / "research/results/a18_microscopy"
OUT = RES / "figs_report"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = ["purewave", "specac", "test_grid_33_11", "test_grid_40_20"]
TITLE = {
    "purewave": "purewave",
    "specac": "specac",
    "test_grid_33_11": "test_grid_33_11",
    "test_grid_40_20": "test_grid_40_20",
}
# значения, стоявшие в fit_lib.GEOMETRY до этой работы
IN_CODE = {
    "purewave": dict(P=25.5, D=10.6),
    "specac": dict(P=24.9, D=14.0),
    "test_grid_33_11": dict(P=31.82, D=11.2),
    "test_grid_40_20": dict(P=38.8, D=18.0),
}

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a85"
GRID = "#e3e3df"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SERIES = [S1, S2, S3, S4]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.edgecolor": GRID,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "text.color": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "figure.dpi": 150,
})


def load(name: str) -> dict:
    return json.loads((RES / name).read_text(encoding="utf-8"))


DATA = {s: load(f"a18_{s}.json") for s in SAMPLES}
CSV = load("a18_csv_measurements.json")
FOCUS = load("a18_test_grid_40_20_focus.json")


def tidy(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis=grid_axis, alpha=0.9)
    ax.set_axisbelow(True)


def save(fig, name: str):
    path = OUT / name
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}")
    return path


# ----------------------------------------------------------------------------------
# Р1/Р2 — сверка источников: период и диаметр
# ----------------------------------------------------------------------------------
def fig_sources(kind: str, fname: str, title: str):
    """Точечная диаграмма: паспорт / обмер владельца / расчёт по снимку / было в коде."""
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ys = np.arange(len(SAMPLES))[::-1]

    for y, s in zip(ys, SAMPLES):
        d = DATA[s]
        passport = d["passport"]["P_um" if kind == "P" else "D_um"]
        code = IN_CODE[s]["P" if kind == "P" else "D"]
        c = CSV[s].get(kind)
        owner = c["all"]["median"] if c else None
        if kind == "P":
            mine, err = d["P_density"]["mean"], d["P_density"]["sem"]
        else:
            mine, err = d["pooled"]["D"]["median"], None

        ax.plot([code], [y], marker="x", ms=8, mew=1.8, color=MUTED, ls="none",
                label="было в коде модели" if y == ys[0] else None, zorder=2)
        ax.plot([passport], [y], marker="s", ms=7, color=S3, ls="none",
                label="паспорт" if y == ys[0] else None, zorder=3)
        if owner is not None:
            ax.plot([owner], [y], marker="D", ms=6.5, color=S2, ls="none",
                    label="обмер владельца (медиана)" if y == ys[0] else None, zorder=4)
        if err:
            ax.errorbar([mine], [y], xerr=[err], fmt="o", ms=8, color=S1,
                        ecolor=S1, elinewidth=1.6, capsize=3,
                        label="расчёт по снимку" if y == ys[0] else None, zorder=5)
        else:
            ax.plot([mine], [y], marker="o", ms=8, color=S1, ls="none",
                    label="расчёт по снимку" if y == ys[0] else None, zorder=5)
        ax.annotate(f"{mine:.2f}", (mine, y), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=7.5, color=INK)
        if abs(code - mine) / mine > 0.03:
            # подпись сбоку от крестика: снизу она налезала бы на ось
            side = 1 if code >= mine else -1
            ax.annotate(f"{code:g}", (code, y), textcoords="offset points",
                        xytext=(9 * side, -1), ha="left" if side > 0 else "right",
                        va="center", fontsize=7.5, color=MUTED)

    ax.set_yticks(ys, [TITLE[s] for s in SAMPLES])
    ax.set_xlabel("мкм")
    ax.set_title(title, loc="left")
    tidy(ax, grid_axis="x")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=4)
    return save(fig, fname)


# ----------------------------------------------------------------------------------
# Р3 — гистограммы периодов: снимок против ручных обмеров владельца
# ----------------------------------------------------------------------------------
def fig_period_hists():
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.0))
    for ax, s in zip(axes.ravel(), SAMPLES):
        d = DATA[s]
        mine = np.asarray(d["pooled"]["periods_um"], float)
        c = CSV[s]["P"]
        owner = np.asarray(c["clean"]["values"] if "values" in c.get("clean", {}) else c["values"], float)
        if c.get("n_outliers"):
            out = set(c["outliers"])
            owner = np.asarray([v for v in c["values"] if v not in out], float)

        lo = min(mine.min(), owner.min())
        hi = max(mine.max(), owner.max())
        bins = np.linspace(lo, hi, 16)
        # плотность, а не счёт: выборки разного объёма (по снимку 9–53, у владельца 47–52)
        ax.hist(owner, bins=bins, density=True, color=S2, alpha=0.45,
                label=f"обмеры владельца (N = {owner.size})", zorder=2)
        ax.hist(mine, bins=bins, density=True, histtype="step", lw=1.8, color=S1,
                label=f"по снимку (N = {mine.size})", zorder=3)

        P = d["P_density"]["mean"]
        ax.axvline(P, color=S1, lw=1.4, ls="--", zorder=4,
                   label=f"средний период {P:.2f} мкм")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.28)
        sm = 100 * d["pooled"]["sigma_P_over_P"]
        so = 100 * np.std(owner, ddof=1) / np.mean(owner)
        ax.set_title(f"{TITLE[s]}\nразброс σ/P: по снимку {sm:.0f} % · у владельца {so:.0f} %",
                     loc="left", fontsize=9)
        ax.set_xlabel("период, мкм")
        ax.set_ylabel("плотность")
        ax.legend(loc="upper left", fontsize=7)
        tidy(ax)
    fig.suptitle("Распределение периода: два независимых измерения", x=0.005, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return save(fig, "r03_period_hist.png")


# ----------------------------------------------------------------------------------
# Р4 — векторы отклонений от идеальной решётки
# ----------------------------------------------------------------------------------
def fig_deviation_vectors():
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 4.8))
    for ax, s in zip(axes.ravel(), SAMPLES):
        d = DATA[s]
        dev = np.asarray(d["pooled"]["dev_lsq_um"], float)
        sd = d["pooled"]["dev"]["std"]
        x = np.arange(1, dev.size + 1)
        ax.axhspan(-sd, sd, color=S1, alpha=0.10, zorder=1)
        ax.axhline(0, color=MUTED, lw=1.0, zorder=2)
        ax.vlines(x, 0, dev, color=S1, lw=1.2, alpha=0.8, zorder=3)
        ax.plot(x, dev, "o", ms=3.4, color=S1, zorder=4)
        P = d["P_density"]["mean"]
        ax.set_title(f"{TITLE[s]}\nσ = {sd:.2f} мкм = {100*sd/P:.0f} % периода",
                     loc="left", fontsize=9)
        ax.set_xlabel("проволока (сквозная нумерация по кадрам)")
        ax.set_ylabel("отклонение, мкм")
        tidy(ax)
    fig.suptitle("Отклонение каждой проволоки от идеальной решётки", x=0.005, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "r04_deviation_vectors.png")


# ----------------------------------------------------------------------------------
# Р5 — гистограммы диаметров
# ----------------------------------------------------------------------------------
def fig_diameter_hists():
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 4.8))
    for ax, s in zip(axes.ravel(), SAMPLES):
        d = DATA[s]
        w = np.asarray(d["pooled"]["widths_um"], float)
        ax.hist(w, bins=14, color=S1, alpha=0.75, zorder=2)
        med = d["pooled"]["D"]["median"]
        pas = d["passport"]["D_um"]
        ax.axvline(med, color=S1, lw=1.6, ls="--", zorder=3)
        ax.axvline(pas, color=S3, lw=1.6, zorder=3)
        c = CSV[s].get("D")
        if c:
            ax.axvline(c["all"]["median"], color=S2, lw=1.6, ls=":", zorder=3)
        ax.set_title(f"{TITLE[s]}\nD = {med:.2f} мкм (N = {d['pooled']['D']['n']}), "
                     f"σ = {d['pooled']['D']['std']:.2f}", loc="left", fontsize=9)
        ax.set_xlabel("ширина проволоки, мкм")
        ax.set_ylabel("число проволок")
        tidy(ax)
    axes[0, 0].plot([], [], color=S1, ls="--", label="медиана по снимку")
    axes[0, 0].plot([], [], color=S3, label="паспорт")
    axes[0, 0].plot([], [], color=S2, ls=":", label="обмер владельца")
    axes[0, 0].legend(loc="upper left")
    fig.suptitle("Распределение диаметра (граница на полувысоте)", x=0.005, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "r05_diameter_hist.png")


# ----------------------------------------------------------------------------------
# Р6 — сводная нерегулярность
# ----------------------------------------------------------------------------------
def fig_irregularity():
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    vals = [100 * DATA[s]["pooled"]["ratio_std_over_P"] for s in SAMPLES]
    ys = np.arange(len(SAMPLES))[::-1]
    ax.barh(ys, vals, height=0.55, color=S1, zorder=2)
    for y, v, s in zip(ys, vals, SAMPLES):
        sd = DATA[s]["pooled"]["dev"]["std"]
        mx = DATA[s]["pooled"]["dev_max_abs_um"]
        ax.annotate(f"{v:.1f} %   (σ = {sd:.2f} мкм, макс. {mx:.1f})",
                    (v, y), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK)
    ax.set_yticks(ys, [TITLE[s] for s in SAMPLES])
    ax.set_xlabel("σ отклонений в процентах периода")
    ax.set_xlim(0, max(vals) * 1.55)
    ax.set_title("Нерегулярность намотки", loc="left")
    tidy(ax, grid_axis="x")
    return save(fig, "r06_irregularity.png")


# ----------------------------------------------------------------------------------
# Р7 — влияние уровня границы
# ----------------------------------------------------------------------------------
def _scan_series(d, key):
    ts = d["threshold_scan"]
    lv, val = [], []
    for k in sorted(ts.keys(), key=float):
        node = ts[k]
        v = node.get(key)
        if isinstance(v, dict):
            v = v.get("mean", v.get("median"))
        if v is not None:
            lv.append(float(k))
            val.append(float(v))
    return np.asarray(lv), np.asarray(val)


def fig_threshold():
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2), sharex=True)
    for i, (key, ax, ttl) in enumerate(
            [("P", axes[0], "Период"), ("D", axes[1], "Диаметр")]):
        for s, col in zip(SAMPLES, SERIES):
            lv, val = _scan_series(DATA[s], key)
            if lv.size == 0:
                continue
            ref = val[np.argmin(np.abs(lv - 0.5))]
            ax.plot(lv, 100 * (val - ref) / ref, "-o", ms=3.5, lw=1.6, color=col,
                    label=TITLE[s])
        ax.axvline(0.5, color=MUTED, lw=1.0, ls="--")
        ax.set_xlabel("уровень границы, доля контраста")
        ax.set_ylabel("отклонение от значения при 0.5, %")
        ax.set_title(ttl, loc="left")
        tidy(ax)
    axes[0].legend(loc="upper left", ncol=2)
    lim = max(abs(np.array(axes[1].get_ylim())).max(), 1)
    axes[0].set_ylim(-lim, lim)
    axes[1].set_ylim(-lim, lim)
    fig.suptitle("Одна и та же шкала: порог двигает диаметр, но не период",
                 x=0.005, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save(fig, "r07_threshold.png")


# ----------------------------------------------------------------------------------
# Р8 — фокус-серия
# ----------------------------------------------------------------------------------
def fig_focus():
    rows = []
    for k, v in FOCUS.items():
        if k.startswith("_"):
            continue
        for fr in v:
            r = fr.get("edge_rise_2080_um")
            if r is None:
                continue
            rows.append((r, fr.get("D_median_um"), fr.get("D20_median_um"),
                         fr.get("P_density_um")))
    rows = [r for r in rows if r[0] is not None]
    rise = np.array([r[0] for r in rows], float)
    summ = FOCUS["_summary"]

    fig, axes = plt.subplots(1, 3, figsize=(7.8, 2.9))
    for ax, (idx, key, ttl, col) in zip(axes, [
            (3, "P_vs_rise", "период", S1),
            (1, "D50_vs_rise", "диаметр, граница 0.5", S3),
            (2, "D20_vs_rise", "диаметр, граница 0.2", S2)]):
        y = np.array([r[idx] if r[idx] is not None else np.nan for r in rows], float)
        m = np.isfinite(y)
        ax.plot(rise[m], y[m], "o", ms=5, color=col, alpha=0.85, zorder=3)
        st = summ[key]
        xs = np.linspace(rise.min(), rise.max(), 20)
        ax.plot(xs, st["intercept"] + st["slope"] * xs, "-", lw=1.8, color=col, zorder=4)
        ax.set_title(f"{ttl}\nr = {st['r']:+.2f}", loc="left")
        ax.set_xlabel("размытие края, мкм")
        ax.set_ylabel("мкм")
        tidy(ax)
    fig.suptitle("Фокус-серия test_grid_40_20: что портится при расфокусировке (19 кадров)",
                 x=0.005, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return save(fig, "r08_focus.png")


# ----------------------------------------------------------------------------------
# Р9 — нуль-модель: беспорядок против дрейфа
# ----------------------------------------------------------------------------------
def fig_null_model():
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.2))
    for ax, key, ttl in [(axes[0], "acf1_periods", "Автокорреляция периодов (лаг 1)"),
                         (axes[1], "sigma_ratio", "σ периодов / σ отклонений")]:
        ys = np.arange(len(SAMPLES))[::-1]
        lo, hi = [], []
        for y, s in zip(ys, SAMPLES):
            nm = DATA[s]["null_model_iid"][key]
            ax.plot([nm["p2_5"], nm["p97_5"]], [y, y], lw=9, color=GRID,
                    solid_capstyle="butt", zorder=2)
            obs = DATA[s]["structure"][f"{key}_mean" if key == "acf1_periods" else "sigma_ratio_mean"]
            ax.plot([obs], [y], "o", ms=8, color=S1, zorder=4)
            ax.annotate(f"{obs:+.2f}" if key == "acf1_periods" else f"{obs:.2f}",
                        (obs, y), xytext=(11, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=8, color=INK, zorder=5)
            lo.append(nm["p2_5"]); hi.append(max(nm["p97_5"], obs))
        pad = 0.22 * (max(hi) - min(lo))
        ax.set_xlim(min(lo) - pad * 0.4, max(hi) + pad * 1.5)
        ax.set_ylim(-0.7, len(SAMPLES) - 0.3)
        ax.set_yticks(ys, [TITLE[s] for s in SAMPLES])
        ax.set_title(ttl, loc="left", fontsize=9.5)
        tidy(ax, grid_axis="x")
    axes[0].plot([], [], lw=9, color=GRID, label="коридор 95 % для случайного джиттера")
    axes[0].plot([], [], "o", ms=8, color=S1, label="измерено")
    axes[0].legend(loc="upper center", bbox_to_anchor=(1.05, -0.18), ncol=2)
    fig.suptitle("Проверка: беспорядок или дрейф намотки", x=0.005, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    return save(fig, "r09_null_model.png")


# ----------------------------------------------------------------------------------
# Р10 — автокорреляционная оценка против метода по границам
# ----------------------------------------------------------------------------------
def fig_acf_vs_border():
    """Точечная диаграмма: 4 серии all-pairs -> палитра blue/orange/aqua/violet
    (проверена `--pairs all`) + разная форма маркера как вторичное кодирование."""
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    cols = [S1, S2, S3, "#4a3aa7"]
    marks = ["o", "s", "^", "D"]
    n_bad = n_tot = 0
    for s, col, mk in zip(SAMPLES, cols, marks):
        d = DATA[s]
        xs, ys = [], []
        for fr, fs in zip(d["frames"], d["frames_summary"]):
            a, m = fr.get("p_autocorr_um"), fs.get("P_density_um")
            if a is None or m is None:
                continue
            xs.append(m)
            ys.append(a)
            n_tot += 1
            if a / m > 1.4 or a / m < 0.7:
                n_bad += 1
        ax.plot(xs, ys, mk, ms=6.5, color=col, alpha=0.9, ls="none",
                mec=SURFACE, mew=0.8, label=TITLE[s], zorder=3)
    lim = [0, 95]
    ax.plot(lim, lim, "-", lw=1.2, color=MUTED, zorder=2)
    ax.plot(lim, [2 * v for v in lim], "--", lw=1.2, color=MUTED, zorder=2)
    ax.annotate("оценки совпадают", (44, 45.5), fontsize=7.5, color=INK2,
                ha="right", va="bottom")
    ax.annotate("завышение вдвое", (30, 61.5), fontsize=7.5, color=INK2,
                ha="left", va="bottom")
    ax.set_xlim(15, 50)
    ax.set_ylim(15, 88)
    ax.set_xlabel("период по границам, мкм")
    ax.set_ylabel("период по автокорреляции, мкм")
    ax.set_title(f"Спектральная оценка промахивается на {n_bad} кадрах из {n_tot}",
                 loc="left")
    tidy(ax, grid_axis="both")
    ax.legend(loc="upper left")
    return save(fig, "r10_acf_vs_border.png")


if __name__ == "__main__":
    print("Рисунки -> research/results/a18_microscopy/figs_report/")
    fig_sources("P", "r01_sources_period.png", "Период: четыре источника")
    fig_sources("D", "r02_sources_diameter.png", "Диаметр: четыре источника")
    fig_period_hists()
    fig_deviation_vectors()
    fig_diameter_hists()
    fig_irregularity()
    fig_threshold()
    fig_focus()
    fig_null_model()
    fig_acf_vs_border()
    print("готово")
