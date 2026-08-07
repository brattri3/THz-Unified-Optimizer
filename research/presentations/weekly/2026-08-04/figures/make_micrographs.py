#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Дифракционные картины двух зон `purewave` — заводской и восстановленной.

⚠ ОСОБЫЙ СТАТУС: как и `make_overlay.py`, работает с исходными данными других зон по
прямому разрешению владельца от 2026-08-04. Никаких измерений по кадрам НЕ производит —
только обрезка и одинаковая для обоих кадров подсветка. Числовые метрики зон уже
посчитаны зоной A (`research/results/two_wgp/a8_diffraction.json`) и цитируются в отчёте.

Что делает: берёт по одному снимку из каждой зоны, вырезает КВАДРАТ вокруг картины
(исходник 4000×1848 — панорамный, по краям пусто; квадрат центрируется на нулевом порядке,
то есть на самом ярком пикселе, — строго геометрический центр кадра сместил бы картину к
краю, потому что линия порядков лежит выше середины), поднимает яркость обоих кадров ОДНИМ
И ТЕМ ЖЕ гамма-преобразованием (иначе на светлом слайде порядки не видны) и кладёт рядом.

⚠ Оговорки, обязательные к переносу в подпись слайда (A8): масштаб кадров НЕ фиксирован —
зум менялся между снимками, поэтому сравнивать по этим двум кадрам расстояния между
порядками нельзя; выборка всего 3 снимка на 3. Кадры иллюстрируют, ЧТО именно измерялось,
а количественный вывод — только из метрик A8, и он имеет статус тенденции.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe research/presentations/weekly/2026-08-04/figures/make_micrographs.py
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
sys.path.insert(0, str(HERE))
from pstyle import INK, INK_2, TITLE  # noqa: E402

SRC = REPO / "research" / "validation" / "micrographs" / "purewave"
GAMMA = 0.55          # одинаково для обоих кадров
CROP = 1150           # сторона вырезаемого квадрата в пикселях исходника
PANELS = [
    ("good 2.jpg", "ЗАВОДСКАЯ зона",
     "порядки — компактные пятна, разделены тёмными промежутками"),
    ("bad 1.jpg", "ВОССТАНОВЛЕННАЯ зона",
     "порядки шире и размытее, между ними — сплошная засветка"),
]


def square_on_zero_order(path, size=900):
    """Квадрат вокруг нулевого порядка + одинаковая для обоих кадров гамма-подсветка.

    Нулевой порядок — самый яркий участок кадра. Ищем его по сильно уменьшенной копии
    (чтобы попасть в центр пятна, а не в случайный «горячий» пиксель матрицы), затем
    вырезаем квадрат `CROP`×`CROP` вокруг найденной точки, прижимая его к границам кадра.
    Строго геометрический центр здесь не годится: исходник панорамный (4000×1848),
    линия порядков лежит выше середины и по краям кадра пусто.
    """
    im = Image.open(path).convert("RGB")
    w, h = im.size

    small = im.convert("L").resize((w // 40, h // 40), Image.BOX)
    g = np.asarray(small, dtype=np.float64)
    iy, ix = np.unravel_index(int(np.argmax(g)), g.shape)
    cx, cy = int((ix + 0.5) * 40), int((iy + 0.5) * 40)

    s = min(CROP, w, h)
    x0 = int(np.clip(cx - s // 2, 0, w - s))
    y0 = int(np.clip(cy - s // 2, 0, h - s))

    im = im.crop((x0, y0, x0 + s, y0 + s)).resize((size, size), Image.LANCZOS)
    a = np.asarray(im, dtype=np.float64) / 255.0
    return np.clip(a ** GAMMA, 0, 1)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.4))
    for ax, (fname, title, sub) in zip(axes, PANELS):
        ax.imshow(square_on_zero_order(SRC / fname))
        ax.set_title(title, color=TITLE, fontsize=13, fontweight="bold", pad=8)
        ax.set_xlabel(sub, color=INK_2, fontsize=9.5, labelpad=6, wrap=True)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor(INK_2)
            sp.set_linewidth(1.0)
    fig.suptitle("Дифракция лазерной указки на решётке: одна указка, два участка образца",
                 color=INK, fontsize=12.5, y=0.98)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.10, wspace=0.06)
    fig.savefig(HERE / "fig6_diffraction_zones.png", dpi=170)
    plt.close(fig)
    print("готово: fig6 (кадры:", ", ".join(p[0] for p in PANELS) + ")")


if __name__ == "__main__":
    main()
