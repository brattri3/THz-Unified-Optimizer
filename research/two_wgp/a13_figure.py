#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A13 — рисунок к отчёту по test_grid_33_11: U(theta) с двумя тенями и закон H5.

Числа НЕ пересчитываются: панель (а) берёт `quality.dynamic_range.per_angle` из
`a13_test_grid_33_11.json`, панель (б) — таблицу семи прежних образцов из
`BRIEF_FOR_P_week_20260728_0804.md` §3.1 плюс новую точку из того же JSON.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "research" / "results" / "two_wgp"
SRC = RES / "a13_test_grid_33_11.json"

# Семь прежних образцов (полный угловой диапазон) — BRIEF_FOR_P §3.1
OLD = [
    ("att-11-16-356", 0.710, 0.373), ("att-11-16-s3", 0.688, 0.371),
    ("att-11-16-s2", 0.688, 0.452),  ("att-11-16-s1", 0.688, 0.528),
    ("specac", 0.562, 0.443),        ("test_grid_40_20", 0.464, 0.586),
    ("purewave (день 1)", 0.416, 0.659),
]


def main():
    d = json.loads(SRC.read_text(encoding="utf-8"))
    per = d["quality"]["dynamic_range"]["per_angle"]
    ang = np.array([p["angle"] for p in per])
    U = np.array([p["U_parseval"] for p in per])
    U_db = 10 * np.log10(U / U.max())

    # Геометрия ИЗМЕРЕНА по микрофото (A14) и отличается от номинальной, на которой
    # строилось предсказание. Показываем и то и другое: старая точка помечена как
    # отозванная, иначе рисунок противоречил бы отчёту.
    ref = json.loads((RES / "a14_refit_measured_geometry.json").read_text(encoding="utf-8"))
    # Диаметр переизмерен с поправкой на расфокусировку (A14b, замечание владельца):
    # 12.55 -> 10.91 мкм. Период и D_eff при этом не меняются (D_phys входит только
    # в знаменатель отношения — проверено повторным фитом).
    # Итоговые числа — из единого пересчёта A16 (поправка A14b отозвана: её
    # экстраполяция переоценивала диаметр). Полная картина по всем образцам —
    # a16_h5_recompute.png; здесь показана только точка этого образца.
    a16 = json.loads((RES / "a16_h5_recompute.json").read_text(encoding="utf-8"))
    r33 = next(r for r in a16["points"] if r["sample"] == "test_grid_33_11")
    dp_new = r33["D_over_P"]
    ratio_new = r33["D_eff_over_D_phys"]
    prereg = 1 - 0.85 * dp_new

    dp_nom = d["D_over_P"]
    ratio_nom = d["fit"]["full_range"]["best"]["D_over_Dphys"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.9))

    # --- (а) угловая кривая ------------------------------------------------
    ax1.plot(ang, U_db, "o-", ms=4, lw=1.2, color="#1f4e79")
    for x, lab, col in ((-90, "тень A", "#c0504d"), (90, "тень B", "#c0504d")):
        ax1.axvline(x, ls="--", lw=1, color=col, alpha=0.55)
    amin = d["quality"]["dynamic_range"]["angle_U_min_parabolic"]
    ax1.axvline(amin, ls=":", lw=1.4, color="#7030a0")
    ax1.annotate(f"минимум {amin:.1f}°", xy=(amin, U_db.min()),
                 xytext=(amin - 78, U_db.min() + 5.5), fontsize=9, color="#7030a0",
                 arrowprops=dict(arrowstyle="->", color="#7030a0", lw=1))
    ax1.text(-90, 1.5, "тень A", ha="center", fontsize=9, color="#c0504d")
    ax1.text(90, 1.5, "тень B", ha="center", fontsize=9, color="#c0504d")
    ax1.set_xlabel("угол оси пропускания θ, град")
    ax1.set_ylabel("U(θ) относительно максимума, дБ")
    ax1.set_title("(а) Пропускание по углу: обе тени сняты\n"
                  f"экстинкция {d['quality']['dynamic_range']['extinction_db']:.1f} дБ, "
                  f"37 углов", fontsize=10.5)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(U_db.min() - 3, 4)

    # --- (б) закон H5 -------------------------------------------------------
    x_old = np.array([o[1] for o in OLD])
    y_old = np.array([o[2] for o in OLD])
    xs = np.linspace(0.33, 0.75, 100)
    ax2.plot(xs, 1 - 0.85 * xs, "-", lw=1.6, color="#404040",
             label="закон H5:  1 − 0.85·D/P")
    ax2.axvspan(0.416, 0.710, color="#dbe5f1", alpha=0.75, zorder=0)
    ax2.text(0.575, 0.735, "диапазон, на котором\nзакон подбирался",
             ha="center", va="top", fontsize=8.5, color="#31538f")
    ax2.plot(x_old, y_old, "o", ms=7, color="#1f4e79", label="7 прежних образцов")
    ax2.plot([dp_nom], [ratio_nom], "o", ms=8, mfc="none", mec="#909090", mew=1.5,
             label=f"номинальная геометрия — ОТОЗВАНО ({ratio_nom:.3f})")
    ax2.annotate("", xy=(dp_new, ratio_new), xytext=(dp_nom, ratio_nom),
                 arrowprops=dict(arrowstyle="->", color="#909090", lw=1.2, ls=":"))
    ax2.plot([dp_new], [prereg], "*", ms=17, mfc="none", mec="#c0504d", mew=1.8,
             label=f"закон H5 при измеренной D/P: {prereg:.3f}")
    ax2.plot([dp_new], [ratio_new], "s", ms=9, color="#c0504d",
             label=f"измерено: {ratio_new:.3f}")
    ax2.annotate("", xy=(dp_new, ratio_new), xytext=(dp_new, prereg),
                 arrowprops=dict(arrowstyle="->", color="#c0504d", lw=1.2))
    ax2.set_xlabel("D/P (плотность решётки)")
    ax2.set_ylabel(r"$D_{eff}/D_{phys}$")
    ax2.set_title("(б) Закон H5 после измерения геометрии по микрофото\n"
                  f"D/P = {dp_new:.3f} (было {dp_nom:.3f}), отклонение "
                  f"{ratio_new - prereg:+.3f}", fontsize=10.5)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8.5, loc="lower left")
    ax2.set_xlim(0.33, 0.75)

    fig.tight_layout()
    out = RES / "a13_test_grid_33_11.png"
    fig.savefig(out, dpi=155)
    print(f"рисунок: {out}")


if __name__ == "__main__":
    main()
