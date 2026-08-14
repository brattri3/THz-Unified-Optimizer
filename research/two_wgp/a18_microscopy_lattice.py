#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A18 — микрофотографии WGP: координаты проволок, период, нерегулярность решётки.

Зона A. Реализует алгоритм владельца (хэндофф ORCH→A от 2026-08-09):

  1. центр проволоки = среднее двух её границ;
  2. проволоки нумеруются по порядку;
  3. по центрам считается средний период;
  4. строится идеальная решётка с этим периодом, привязанная к ПЕРВОЙ проволоке;
  5. вектор отклонений реальных центров от идеальных = мера нерегулярности.

Дополнительно (не вместо, а рядом с методом владельца):
  * привязка решётки по МНК — как контроль к привязке «за первую проволоку»;
  * скан порога границы f = 0.2…0.8 — прямая проверка оговорки ORCH: сдвиг границы
    обязан двигать ДИАМЕТР и не двигать ПЕРИОД;
  * структура отклонений: тренд, автокорреляция, тест ACF(периодов) = -0.5,
    Монте-Карло нуль-модель независимого джиттера.

Подкоманды
----------
    csv      разбор CSV-обмеров владельца + робастная статистика + вопрос «+0.2 мкм»
    frame    один кадр, подробная диагностика (наклон, профиль, порог)
    sample   все основные кадры одного образца -> JSON + PNG
    all      четыре образца (основные кадры) -> JSON + PNG
    focus    фокус-серии test_grid_40_20 (запускать последней)

Запуск ТОЛЬКО из корня репозитория:
    .venv\\Scripts\\python.exe research/two_wgp/a18_microscopy_lattice.py all
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

OUT = REPO / "research" / "results" / "a18_microscopy"

# --------------------------------------------------------------------------------------
# Константы поставки (data_pool/MICROSCOPY.md — единственный источник калибровки)
# --------------------------------------------------------------------------------------
CAL_UM_PER_PX = 0.126739

#: набор -> паспорт (D, P) в мкм, префикс CSV, режим (ожидаемая полярность),
#: список основных кадров (относительно data_pool/<набор>/microscopy/raw)
SAMPLES = {
    "purewave": dict(
        csv_prefix="purewave-10-25",
        D_pass=10.0, P_pass=25.0,
        expect_wire="bright",          # негатив применён
        raw_glob="purewave10_8_2026_*.bmp",
    ),
    "specac": dict(
        csv_prefix="specac-10-25",
        D_pass=10.0, P_pass=25.0,
        expect_wire="bright",
        raw_glob="Image10_8_2026_*.bmp",
    ),
    "test_grid_33_11": dict(
        csv_prefix="wgp-11-33",
        D_pass=11.0, P_pass=33.0,
        expect_wire="bright",
        raw_glob="Image7_8_2026_*.bmp",
    ),
    "test_grid_40_20": dict(
        csv_prefix="wgp-20-40",
        D_pass=20.0, P_pass=40.0,
        expect_wire="dark",            # негатива нет
        raw_glob="multifocus/multifocus_*.bmp",
    ),
}

#: прежние значения research/experiments/fit_lib.py:GEOMETRY (сверка, п.6 ТЗ)
GEOMETRY_OLD = {
    "purewave":        dict(P_um=25.5, D_phys_um=10.6),
    "specac":          dict(P_um=24.9, D_phys_um=14.0),
    "test_grid_33_11": dict(P_um=31.82, D_phys_um=11.2),   # P из A13/A16, D — решение владельца
    "test_grid_40_20": dict(P_um=38.8, D_phys_um=18.0),
}

#: прежние FFT-оценки периода (зона A, A14/A16), мкм — для п.8 ТЗ
FFT_OLD_P_UM = {
    "purewave": [24.72, 26.10],        # неоднозначность 5.6 %, HANDOFF_A_20260807 §5.1
    "specac": [24.91],
    "test_grid_33_11": [31.82],
    "test_grid_40_20": [38.87],
}


# ======================================================================================
# CSV — разбор выгрузки ПО микроскопа (cp1251, ';', десятичная запятая, \xa0 в заголовках)
# ======================================================================================
def _norm_header(s: str) -> str:
    return s.replace("\xa0", " ").replace(" ", " ").strip().lower()


def read_measure_csv(path: Path):
    """Читает `_D.csv` / `_P.csv`: возвращает список значений колонки «Длина, мкм»."""
    raw = path.read_bytes().decode("cp1251")
    lines = [ln for ln in raw.splitlines() if ln.strip()]
    head = [_norm_header(c) for c in lines[0].split(";")]
    try:
        col = next(i for i, c in enumerate(head) if c.startswith("длина"))
    except StopIteration:
        raise ValueError(f"{path}: колонка «Длина» не найдена, заголовки={head}")
    vals = []
    for ln in lines[1:]:
        parts = ln.split(";")
        if len(parts) <= col:
            continue
        tok = parts[col].strip()
        if tok in ("", "-"):
            continue
        vals.append(float(tok.replace(",", ".")))
    return np.asarray(vals, dtype=float)


def read_pstat_csv(path: Path):
    """Читает `_P_stat.csv` — гистограмму ПО микроскопа. Не независимый источник."""
    raw = path.read_bytes().decode("cp1251")
    lines = [ln for ln in raw.splitlines() if ln.strip()]
    head = [_norm_header(c) for c in lines[0].split(";")]

    def idx(prefix):
        for i, c in enumerate(head):
            if c.startswith(prefix):
                return i
        return None

    # заголовки: '#','Название','Минимум X, мкм','Максимум X, мкм','Количество',...
    i_lo, i_hi, i_n = idx("минимум x"), idx("максимум x"), idx("количество")
    if i_lo is None or i_hi is None or i_n is None:
        i_lo, i_hi, i_n = 2, 3, 4
    rows = []
    for ln in lines[1:]:
        p = ln.split(";")
        if len(p) <= max(i_lo, i_hi, i_n):
            continue
        try:
            rows.append(dict(lo=float(p[i_lo].replace(",", ".")),
                             hi=float(p[i_hi].replace(",", ".")),
                             n=int(float(p[i_n].replace(",", ".")))))
        except ValueError:
            continue
    return rows


def robust_stats(x, name=""):
    """Классические + робастные оценки. MAD нормирован к σ (×1.4826)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(name=name, n=0)
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    q1, q3 = np.percentile(x, [25, 75])
    return dict(
        name=name, n=int(x.size),
        mean=float(np.mean(x)), std=float(np.std(x, ddof=1)) if x.size > 1 else 0.0,
        median=med, mad=mad, mad_sigma=1.4826 * mad,
        q1=float(q1), q3=float(q3), iqr=float(q3 - q1),
        min=float(x.min()), max=float(x.max()), ptp=float(np.ptp(x)),
    )


def mad_outliers(x, k=3.5):
    """Маска выбросов по модифицированному z-счёту Иглвича — Хоглина."""
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.zeros_like(x, dtype=bool), np.zeros_like(x)
    z = 0.6745 * (x - med) / mad
    return np.abs(z) > k, z


# ======================================================================================
# Изображение: полярность, наклон, профиль
# ======================================================================================
def load_gray(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("L"), dtype=np.float64)


def detect_polarity(gray: np.ndarray):
    """Определяет, какая фаза — проволока, ПО ФАКТУ, а не по таблице режимов.

    Критерий: промежуток — открытая апертура на просвет, поле там ровное (часто
    клиппировано), проволока круглая и бликует -> внутри неё разброс яркости выше.
    Возвращает (+1 если проволока светлая, -1 если тёмная, диагностика).
    """
    lo, hi = np.percentile(gray, [2, 98])
    thr = 0.5 * (lo + hi)
    bright = gray > thr
    dark = ~bright
    # «ядра» фаз: убираем переходные пиксели (в пределах 20 % контраста от порога)
    band = 0.20 * (hi - lo)
    core_b = gray > thr + band
    core_d = gray < thr - band
    sb = float(gray[core_b].std()) if core_b.sum() > 100 else np.inf
    sd = float(gray[core_d].std()) if core_d.sum() > 100 else np.inf
    pol = 1 if sb > sd else -1
    return pol, dict(std_core_bright=sb, std_core_dark=sd,
                     frac_bright=float(bright.mean()), frac_dark=float(dark.mean()),
                     p2=float(lo), p98=float(hi),
                     wire=("bright" if pol > 0 else "dark"))


def _xcorr_shift(a, b, max_lag):
    """Сдвиг b относительно a (в px, субпиксельно) по кросс-корреляции."""
    a = a - a.mean()
    b = b - b.mean()
    n = a.size
    c = np.correlate(b, a, mode="full")           # индекс n-1 == нулевой лаг
    lags = np.arange(-(n - 1), n)
    m = np.abs(lags) <= max_lag
    c, lags = c[m], lags[m]
    k = int(np.argmax(c))
    if 0 < k < c.size - 1:
        y0, y1, y2 = c[k - 1], c[k], c[k + 1]
        den = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    else:
        delta = 0.0
    return float(lags[k] + delta)


def sheared_profile(gray, slope, row_step=1):
    """Профиль яркости, усреднённый ВДОЛЬ проволок с учётом их наклона.

    slope — сдвиг x на единицу y (px/px). Возвращает (profile, x0, valid_slice).
    """
    H, W = gray.shape
    y0 = (H - 1) / 2.0
    xs = np.arange(W, dtype=np.float64)
    acc = np.zeros(W)
    cnt = np.zeros(W)
    for yy in range(0, H, row_step):
        sh = slope * (yy - y0)
        src = xs + sh
        val = np.interp(src, xs, gray[yy], left=np.nan, right=np.nan)
        m = np.isfinite(val)
        acc[m] += val[m]
        cnt[m] += 1
    full = cnt.max()
    valid = cnt >= full            # только колонки с полным покрытием
    prof = np.full(W, np.nan)
    prof[valid] = acc[valid] / cnt[valid]
    i0 = int(np.argmax(valid))
    i1 = W - int(np.argmax(valid[::-1]))
    return prof[i0:i1], i0


def estimate_slope(gray, nband=8, row_step=2):
    """Наклон проволок: грубо — по сдвигу профилей горизонтальных полос,
    точно — максимизацией контраста сдвинутого профиля."""
    H, W = gray.shape
    edges = np.linspace(0, H, nband + 1).astype(int)
    ys, profs = [], []
    for k in range(nband):
        b = gray[edges[k]:edges[k + 1]]
        profs.append(b.mean(axis=0))
        ys.append(0.5 * (edges[k] + edges[k + 1]))
    ys = np.asarray(ys)
    max_lag = max(4, W // 10)
    sh = np.array([_xcorr_shift(profs[0], p, max_lag) for p in profs])
    A = np.vstack([ys - ys.mean(), np.ones_like(ys)]).T
    coarse = float(np.linalg.lstsq(A, sh, rcond=None)[0][0])

    # точная настройка: максимум СКО профиля (резкость полос)
    span = max(0.004, 0.5 * abs(coarse))
    grid = np.linspace(coarse - span, coarse + span, 41)
    best, best_s = -np.inf, coarse
    for s in grid:
        p, _ = sheared_profile(gray, s, row_step=row_step)
        v = float(np.nanstd(p))
        if v > best:
            best, best_s = v, s
    # второй, более мелкий проход
    step = grid[1] - grid[0]
    grid2 = np.linspace(best_s - step, best_s + step, 21)
    for s in grid2:
        p, _ = sheared_profile(gray, s, row_step=row_step)
        v = float(np.nanstd(p))
        if v > best:
            best, best_s = v, s
    return float(best_s), dict(coarse=coarse, contrast=best,
                               tilt_deg=float(np.degrees(np.arctan(best_s))))


# ======================================================================================
# Границы проволок и центры (алгоритм владельца, шаги 1–2)
# ======================================================================================
def _smooth(p, sigma):
    if sigma <= 0:
        return p.copy()
    r = int(np.ceil(3 * sigma))
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    pad = np.pad(p, r, mode="edge")
    return np.convolve(pad, k, mode="valid")


def _rolling_env(p, win):
    """Огибающие min/max скользящим окном (борьба с виньетированием кадра)."""
    from numpy.lib.stride_tricks import sliding_window_view
    win = int(max(3, win))
    half = win // 2
    pad = np.pad(p, half, mode="edge")
    w = sliding_window_view(pad, 2 * half + 1)[:p.size]
    return _smooth(w.min(axis=1), win / 4.0), _smooth(w.max(axis=1), win / 4.0)


def _segments(mask, min_len):
    """Список (начало, конец) непрерывных участков True длиной >= min_len."""
    m = np.asarray(mask, dtype=np.int8)
    d = np.diff(np.concatenate(([0], m, [0])))
    starts = np.nonzero(d == 1)[0]
    ends = np.nonzero(d == -1)[0]            # полуинтервал [start, end)
    keep = (ends - starts) >= min_len
    return list(zip(starts[keep], ends[keep]))


def _cross(p, thr, i_from, i_to):
    """Первое пересечение уровня thr при движении от i_from к i_to (субпиксельно)."""
    step = 1 if i_to > i_from else -1
    idx = np.arange(i_from, i_to, step)
    if idx.size < 2:
        return None
    v = p[idx] - thr
    s = np.sign(v)
    s[s == 0] = s[0] if s[0] != 0 else 1
    ch = np.nonzero(np.diff(s) != 0)[0]
    if ch.size == 0:
        return None
    k = ch[0]
    a, b = idx[k], idx[k + 1]
    va, vb = p[a] - thr, p[b] - thr
    if va == vb:
        return float(a)
    return float(a + (b - a) * va / (va - vb))


def autocorr_period(p, pmin_px, pmax_px):
    """Оценка периода по автокорреляции профиля (только для выбора окна огибающей)."""
    q = p - p.mean()
    ac = np.correlate(q, q, mode="full")[q.size - 1:]
    ac /= ac[0]
    lo, hi = int(pmin_px), min(int(pmax_px), ac.size - 2)
    if hi <= lo + 2:
        return np.nan
    seg = ac[lo:hi]
    k = int(np.argmax(seg)) + lo
    if 0 < k < ac.size - 1:
        y0, y1, y2 = ac[k - 1], ac[k], ac[k + 1]
        den = (y0 - 2 * y1 + y2)
        d = (0.5 * (y0 - y2) / den) if den != 0 else 0.0
        # плоская вершина -> знаменатель ~0 и поправка улетает; ограничиваем полпикселем
        k = k + float(np.clip(d, -0.5, 0.5))
    return float(k)


def orient_profile(prof, pol, p_approx_px, smooth_px=2.0):
    """Ориентирует профиль (проволока = максимум) и считает огибающие один раз."""
    p = _smooth(prof * pol, smooth_px)
    lo, hi = _rolling_env(p, win=1.25 * p_approx_px)
    return p, lo, hi


def find_wires(p, lo, hi, p_approx_px, level=0.5, gap_level=0.15, min_frac=0.06):
    """Границы и центры проволок — двухстадийно, устойчиво к бликам ВНУТРИ проволоки.

    Стадия 1: находим «ядра промежутка» — участки, где профиль глубоко внизу
      (`p < lo + gap_level*(hi-lo)`). Промежуток — открытая апертура, он всегда
      чистый; проволока может иметь внутри бликовую полосу (`test_grid_40_20`),
      от которой одностадийный порог рассыпается.
    Стадия 2: между двумя соседними ядрами промежутка лежит ровно одна проволока.
      Её границы ищем как ПЕРВОЕ пересечение уровня
      `gap + level*(wire - gap)` при движении ОТ промежутка ВНУТРЬ проволоки —
      то есть по монотонному склону, куда внутренняя структура не достаёт.
    Центр = среднее двух границ (шаг 1 алгоритма владельца).
    """
    contrast = hi - lo
    gap_mask = p < lo + gap_level * contrast
    segs = _segments(gap_mask, min_len=max(3, int(0.08 * p_approx_px)))
    empty = dict(centers=np.array([]), widths=np.array([]),
                 left=np.array([]), right=np.array([]), n=0)
    if len(segs) < 2:
        return empty

    left, right, wlev, cmass = [], [], [], []
    for (a0, a1), (b0, b1) in zip(segs[:-1], segs[1:]):
        # Стартуем поиск границы из ГЛУБИНЫ промежутка (его минимума), а не с края
        # сегмента: на краю профиль уже может оказаться выше низкого порога, и тогда
        # пересечение не находится — проволока теряется. Проверено на
        # test_grid_40_20/multifocus_3: при f=0.50 терялась одна проволока из семи.
        i_a = a0 + int(np.argmin(p[a0:a1]))
        i_b = b0 + int(np.argmin(p[b0:b1]))
        if i_b - i_a < max(3, int(0.05 * p_approx_px)):
            continue
        inner = p[i_a:i_b + 1]
        core = i_a + int(np.argmax(inner))
        wire_level = float(np.percentile(inner, 95))
        gl = 0.5 * (float(np.median(p[a0:a1])) + float(np.median(p[b0:b1])))
        thr = gl + level * (wire_level - gl)
        xl = _cross(p, thr, i_a, core + 1)
        xr = _cross(p, thr, i_b, core - 1)
        if xl is None or xr is None or xr <= xl:
            continue
        # контрольная оценка центра: центр тяжести «затемнения» между промежутками
        wgt = np.clip(inner - gl, 0, None)
        cm = (float(np.sum(wgt * np.arange(i_a, i_b + 1)) / np.sum(wgt))
              if np.sum(wgt) > 0 else 0.5 * (xl + xr))
        left.append(xl)
        right.append(xr)
        wlev.append(wire_level - gl)
        cmass.append(cm)

    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    wlev = np.asarray(wlev, dtype=float)
    cmass = np.asarray(cmass, dtype=float)
    if left.size == 0:
        return empty
    widths = right - left
    keep = widths > min_frac * p_approx_px
    left, right, widths, wlev, cmass = (left[keep], right[keep], widths[keep],
                                        wlev[keep], cmass[keep])
    if left.size == 0:
        return empty
    # Слишком узкие находки — не проволоки: это обрезанная проволока у края кадра либо
    # слабый провал яркости внутри расфокусированной проволоки. Проволока этой решётки
    # имеет вполне определённую ширину, и половины от медианной у неё быть не может.
    # Такие находки УДАЛЯЮТСЯ (не просто помечаются), иначе они входят в решётку как
    # лишние узлы и занижают период.
    med_w0 = float(np.median(widths))
    narrow = widths < 0.55 * med_w0
    n_narrow_inner = int(narrow[1:-1].sum()) if narrow.size > 2 else 0
    if narrow.any():
        left, right, widths, wlev, cmass = (left[~narrow], right[~narrow], widths[~narrow],
                                            wlev[~narrow], cmass[~narrow])
    if left.size == 0:
        return empty
    centers = 0.5 * (left + right)
    # признак ненадёжной проволоки: слабый контраст к промежутку (расфокусировка)
    # или слипшаяся пара (ширина ~2 D)
    med_c = float(np.median(wlev)) if wlev.size else 0.0
    med_w = float(np.median(widths))
    ok = ((wlev > 0.35 * med_c) & (widths < 1.6 * med_w) & (widths < 0.95 * p_approx_px))
    return dict(centers=centers, widths=widths, left=left, right=right,
                centroids=cmass, contrast_wire=wlev, ok=ok,
                n=int(centers.size), n_bad=int((~ok).sum()),
                n_narrow_dropped=int(narrow.sum()), n_narrow_inner=n_narrow_inner)


# ======================================================================================
# Решётка и нерегулярность (шаги 3–5 алгоритма владельца)
# ======================================================================================
def lattice_analysis(centers_um):
    """Средний период, идеальная решётка, вектор отклонений + структура отклонений."""
    x = np.asarray(centers_um, dtype=float)
    n = x.size
    res = dict(n_wires=int(n))
    if n < 3:
        return res
    idx = np.arange(n, dtype=float)
    d = np.diff(x)
    res["periods"] = d.tolist()
    res["period_stats"] = robust_stats(d, "period_um")

    # --- шаг 3-5 владельца: средний период + привязка за ПЕРВУЮ проволоку -------------
    P_mean = float(d.mean())
    ideal_owner = x[0] + idx * P_mean
    dev_owner = x - ideal_owner
    # --- контроль: МНК-решётка (свободные период и сдвиг) -----------------------------
    A = np.vstack([idx, np.ones(n)]).T
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    P_lsq, x0_lsq = float(coef[0]), float(coef[1])
    dev_lsq = x - (x0_lsq + idx * P_lsq)

    res["P_mean_um"] = P_mean
    res["P_median_um"] = float(np.median(d))
    res["P_lsq_um"] = P_lsq
    res["dev_owner_um"] = dev_owner.tolist()
    res["dev_lsq_um"] = dev_lsq.tolist()
    res["dev_owner"] = dict(
        std=float(np.std(dev_owner, ddof=1)), max_abs=float(np.abs(dev_owner).max()),
        rms=float(np.sqrt(np.mean(dev_owner ** 2))),
        ratio_std_over_P=float(np.std(dev_owner, ddof=1) / P_mean),
    )
    res["dev_lsq"] = dict(
        std=float(np.std(dev_lsq, ddof=1)), max_abs=float(np.abs(dev_lsq).max()),
        rms=float(np.sqrt(np.mean(dev_lsq ** 2))),
        ratio_std_over_P=float(np.std(dev_lsq, ddof=1) / P_lsq),
    )

    # --- структура отклонений ---------------------------------------------------------
    # линейный тренд отклонений владельца (равносилен смещению среднего периода)
    A2 = np.vstack([idx, np.ones(n)]).T
    b, *_ = np.linalg.lstsq(A2, dev_owner, rcond=None)
    fit = A2 @ b
    resid = dev_owner - fit
    dof = n - 2
    s2 = float(resid @ resid / dof) if dof > 0 else np.nan
    cov = s2 * np.linalg.inv(A2.T @ A2) if dof > 0 else np.full((2, 2), np.nan)
    slope, slope_se = float(b[0]), float(np.sqrt(cov[0, 0]))
    res["dev_trend"] = dict(slope_um_per_wire=slope, slope_se=slope_se,
                            t=float(slope / slope_se) if slope_se > 0 else np.nan,
                            dof=int(dof))

    def acf1(v):
        v = np.asarray(v, float)
        v = v - v.mean()
        if v.size < 3 or np.all(v == 0):
            return np.nan
        return float((v[:-1] @ v[1:]) / (v @ v))

    res["acf1_dev_lsq"] = acf1(dev_lsq)
    res["acf1_periods"] = acf1(d)
    res["sigma_P_over_sigma_dev"] = (
        float(np.std(d, ddof=1) / np.std(dev_lsq, ddof=1))
        if np.std(dev_lsq, ddof=1) > 0 else np.nan)
    return res


def null_model_iid(n, nmc=4000, seed=0):
    """Нуль-модель: центры = идеальная решётка + независимый гауссов джиттер σ=1.

    Даёт распределение ACF(периодов) и σ_P/σ_dev при том же N — иначе малое N
    само по себе смещает эти оценки, и сравнивать с «-0.5» и «√2» нельзя.
    """
    rng = np.random.default_rng(seed)
    idx = np.arange(n, dtype=float)
    A = np.vstack([idx, np.ones(n)]).T
    pinv = np.linalg.pinv(A)
    acfs, ratios = [], []
    for _ in range(nmc):
        x = idx * 100.0 + rng.standard_normal(n)
        d = np.diff(x)
        dev = x - A @ (pinv @ x)
        v = d - d.mean()
        acfs.append((v[:-1] @ v[1:]) / (v @ v))
        sd = np.std(dev, ddof=1)
        ratios.append(np.std(d, ddof=1) / sd if sd > 0 else np.nan)
    acfs = np.asarray(acfs)
    ratios = np.asarray(ratios, dtype=float)
    return dict(n=n,
                acf1_periods=dict(mean=float(acfs.mean()),
                                  p2_5=float(np.percentile(acfs, 2.5)),
                                  p97_5=float(np.percentile(acfs, 97.5))),
                sigma_ratio=dict(mean=float(np.nanmean(ratios)),
                                 p2_5=float(np.nanpercentile(ratios, 2.5)),
                                 p97_5=float(np.nanpercentile(ratios, 97.5))))


# ======================================================================================
# Обработка кадра
# ======================================================================================
LEVELS = (0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)


def process_frame(path: Path, p_pass_um: float, cal=CAL_UM_PER_PX,
                  levels=LEVELS, row_step=1, verbose=False):
    gray = load_gray(path)
    pol, poldiag = detect_polarity(gray)
    slope, slopediag = estimate_slope(gray, row_step=2)
    prof, x_off = sheared_profile(gray, slope, row_step=row_step)
    prof = prof[np.isfinite(prof)]

    p_pass_px = p_pass_um / cal
    p_ac = autocorr_period(prof * pol, 0.4 * p_pass_px, 2.0 * p_pass_px)
    # автокорреляция используется ТОЛЬКО как масштаб окна огибающей; если она
    # уехала от паспорта больше чем в 1.5 раза — берём паспорт, чтобы окно не
    # развалило сегментацию (само измерение периода от этого числа не зависит).
    if np.isfinite(p_ac) and (1 / 1.5) <= p_ac / p_pass_px <= 1.5:
        p_approx = p_ac
    else:
        p_approx = p_pass_px

    out = dict(
        file=str(path.relative_to(REPO)).replace("\\", "/"),
        width_px=int(gray.shape[1]), height_px=int(gray.shape[0]),
        polarity=poldiag, slope=slopediag, slope_px_per_px=slope,
        p_autocorr_px=float(p_ac) if np.isfinite(p_ac) else None,
        p_autocorr_um=float(p_ac * cal) if np.isfinite(p_ac) else None,
        p_passport_px=float(p_pass_px),
        levels={},
    )
    p_or, lo_env, hi_env = orient_profile(prof, pol, p_approx)
    for lv in levels:
        w = find_wires(p_or, lo_env, hi_env, p_approx, level=lv)
        if w["n"] < 3:
            out["levels"][f"{lv:.2f}"] = dict(n_wires=int(w["n"]))
            continue
        c_um = w["centers"] * cal
        d_um = w["widths"] * cal
        rec = lattice_analysis(c_um)
        rec["centers_px"] = w["centers"].tolist()
        rec["centers_um"] = c_um.tolist()
        rec["widths_um"] = d_um.tolist()
        rec["width_stats"] = robust_stats(d_um, "D_um")
        rec["n_bad_wires"] = int(w["n_bad"])
        rec["n_narrow_dropped"] = int(w.get("n_narrow_dropped", 0))
        rec["n_narrow_inner"] = int(w.get("n_narrow_inner", 0))
        rec["left_um"] = (w["left"] * cal).tolist()
        rec["right_um"] = (w["right"] * cal).tolist()
        # --- средний период по ПЛОТНОСТИ проволок ------------------------------------
        # Устраняет смещение от слипшихся проволок: там, где нерегулярность сводит две
        # проволоки вплотную, детектор видит одну «двойную». Отбрасывать такой кадр
        # нельзя — это выбрасывает именно короткие периоды и завышает средний.
        # Считаем кратность каждой «проволоки» по её ширине и делим пролёт на число
        # межпроволочных промежутков.
        med_w = float(np.median(d_um[w["ok"]])) if w["ok"].any() else float(np.median(d_um))
        units = np.maximum(1, np.round(d_um / med_w)).astype(int)
        n_units = int(units.sum())
        x_first = float(w["left"][0] * cal) + 0.5 * med_w
        x_last = float(w["right"][-1] * cal) - 0.5 * med_w
        rec["n_units"] = n_units
        rec["n_merged"] = int((units > 1).sum())
        rec["P_density_um"] = ((x_last - x_first) / (n_units - 1)) if n_units > 1 else None
        rec["ok_mask"] = w["ok"].astype(int).tolist()
        rec["contrast_wire"] = np.asarray(w["contrast_wire"]).tolist()
        # контроль: та же решётка, но центр = центр тяжести затемнения
        cm_um = w["centroids"] * cal
        lc = lattice_analysis(cm_um)
        rec["centroid_control"] = dict(
            P_mean_um=lc.get("P_mean_um"),
            dev_std_um=(lc.get("dev_lsq", {}) or {}).get("std"),
            period_std_um=(lc.get("period_stats", {}) or {}).get("std"),
            centers_um=cm_um.tolist())
        # диаметр только по надёжным проволокам
        if w["ok"].any():
            rec["width_stats_ok"] = robust_stats(d_um[w["ok"]], "D_um_ok")
        rec["duty_D_over_P"] = (float(np.median(d_um) / rec["P_mean_um"])
                                if rec.get("P_mean_um") else None)
        out["levels"][f"{lv:.2f}"] = rec
    # профиль для рисунка (проредить)
    out["_prof"] = prof
    out["_pol"] = pol
    out["_p_approx"] = p_approx
    if verbose:
        base = out["levels"].get("0.50", {})
        print(f"  {path.name}: wire={poldiag['wire']} tilt={slopediag['tilt_deg']:+.3f}deg "
              f"N={base.get('n_wires', 0)} P={base.get('P_mean_um', float('nan')):.3f} "
              f"D={base.get('width_stats', {}).get('median', float('nan')):.3f}")
    return out


# ======================================================================================
# Рисунки
# ======================================================================================
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_frame(res, out_png: Path, cal=CAL_UM_PER_PX, title=""):
    plt = _mpl()
    prof = res["_prof"]
    pol = res["_pol"]
    lv = res["levels"].get("0.50", {})
    fig, ax = plt.subplots(3, 1, figsize=(11, 9))

    x_um = np.arange(prof.size) * cal
    ax[0].plot(x_um, prof * pol, lw=0.8, color="0.3")
    if lv.get("n_wires", 0) >= 3:
        c = np.asarray(lv["centers_um"])
        for i, cc in enumerate(c):
            ax[0].axvline(cc, color="tab:red", lw=0.8, alpha=0.8)
        ax[0].set_title(f"{title}\nпрофиль (проволока = максимум), центры f=0.50, "
                        f"N={lv['n_wires']}, наклон {res['slope']['tilt_deg']:+.3f}°")
    ax[0].set_xlabel("x, мкм")
    ax[0].set_ylabel("яркость (ориент.)")

    if lv.get("n_wires", 0) >= 3:
        dev = np.asarray(lv["dev_owner_um"])
        devl = np.asarray(lv["dev_lsq_um"])
        k = np.arange(dev.size)
        ax[1].plot(k, dev * 1000, "o-", label="привязка за 1-ю (владелец)")
        ax[1].plot(k, devl * 1000, "s--", label="МНК-решётка")
        ax[1].axhline(0, color="0.6", lw=0.7)
        ax[1].set_xlabel("номер проволоки")
        ax[1].set_ylabel("отклонение, нм")
        ax[1].legend(fontsize=8)
        ax[1].set_title(f"вектор отклонений: СКО(МНК) = {lv['dev_lsq']['std']*1000:.0f} нм, "
                        f"СКО/P = {lv['dev_lsq']['ratio_std_over_P']*100:.2f} %")

        # зависимость P и D от порога
        lvs, Ps, Ds = [], [], []
        for k_, v in sorted(res["levels"].items()):
            if v.get("n_wires", 0) >= 3:
                lvs.append(float(k_))
                Ps.append(v["P_mean_um"])
                Ds.append(v["width_stats"]["median"])
        ax2 = ax[2]
        ax2.plot(lvs, Ps, "o-", color="tab:blue", label="средний период P")
        ax2.set_xlabel("уровень порога границы f")
        ax2.set_ylabel("P, мкм", color="tab:blue")
        ax2b = ax2.twinx()
        ax2b.plot(lvs, Ds, "s-", color="tab:orange", label="медианный диаметр D")
        ax2b.set_ylabel("D, мкм", color="tab:orange")
        ax2.set_title("проверка оговорки ORCH: порог двигает D, не двигает P")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def plot_sample(agg, out_png: Path, title=""):
    plt = _mpl()
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    per = np.asarray(agg["pooled"]["periods_um"])
    dev = np.asarray(agg["pooled"]["dev_lsq_um"])
    wid = np.asarray(agg["pooled"]["widths_um"])

    ax[0, 0].hist(per, bins=max(8, int(np.sqrt(per.size)) * 2), color="tab:blue", alpha=.8)
    ax[0, 0].axvline(np.median(per), color="k", lw=1.2)
    ax[0, 0].set_title(f"периоды по снимкам, N={per.size}\n"
                       f"мед={np.median(per):.3f} мкм, MADσ="
                       f"{1.4826*np.median(np.abs(per-np.median(per))):.3f}")
    ax[0, 0].set_xlabel("P, мкм")

    ax[0, 1].hist(dev * 1000, bins=max(8, int(np.sqrt(max(dev.size, 1))) * 2),
                  color="tab:green", alpha=.8)
    ax[0, 1].set_title(f"отклонения от МНК-решётки, N={dev.size}\n"
                       f"СКО={np.std(dev, ddof=1)*1000:.0f} нм")
    ax[0, 1].set_xlabel("отклонение, нм")

    if wid.size:
        ax[1, 0].hist(wid, bins=max(8, int(np.sqrt(wid.size)) * 2), color="tab:orange", alpha=.8)
        ax[1, 0].axvline(np.median(wid), color="k", lw=1.2)
        ax[1, 0].set_title(f"диаметры (f=0.50), N={wid.size}\nмед={np.median(wid):.3f} мкм")
        ax[1, 0].set_xlabel("D, мкм")

    # период кадр за кадром
    fr = agg["frames_summary"]
    xs = np.arange(len(fr))
    ys = [f["P_mean_um"] for f in fr]
    es = [f["P_sem_um"] for f in fr]
    ax[1, 1].errorbar(xs, ys, yerr=es, fmt="o-")
    ax[1, 1].set_xticks(xs)
    ax[1, 1].set_xticklabels([Path(f["file"]).stem[-12:] for f in fr],
                             rotation=45, ha="right", fontsize=7)
    ax[1, 1].set_ylabel("P, мкм")
    ax[1, 1].set_title("средний период по кадрам")
    fig.suptitle(title)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


# ======================================================================================
# Подкоманды
# ======================================================================================
def cmd_csv(args):
    OUT.mkdir(parents=True, exist_ok=True)
    result = {}
    for name, cfg in SAMPLES.items():
        base = REPO / "data_pool" / name / "microscopy" / "measurements"
        rec = dict(passport=dict(D_um=cfg["D_pass"], P_um=cfg["P_pass"]))
        for kind in ("D", "P"):
            f = base / f"{cfg['csv_prefix']}_{kind}.csv"
            if not f.exists():
                rec[kind] = None
                continue
            v = read_measure_csv(f)
            out, z = mad_outliers(v)
            rec[kind] = dict(
                file=str(f.relative_to(REPO)).replace("\\", "/"),
                values=v.tolist(),
                all=robust_stats(v, kind),
                clean=robust_stats(v[~out], kind + "_clean"),
                n_outliers=int(out.sum()),
                outliers=v[out].tolist(),
            )
        fstat = base / f"{cfg['csv_prefix']}_P_stat.csv"
        if fstat.exists():
            rows = read_pstat_csv(fstat)
            rec["P_stat"] = dict(bins=rows, sum_n=int(sum(r["n"] for r in rows)))
        # вопрос «+0.2 мкм»
        if rec.get("D"):
            for tag in ("all", "clean"):
                s = rec["D"][tag]
                rec.setdefault("D_vs_passport", {})[tag] = dict(
                    delta_mean=s["mean"] - cfg["D_pass"],
                    delta_median=s["median"] - cfg["D_pass"],
                    ratio_mean=s["mean"] / cfg["D_pass"],
                    ratio_median=s["median"] / cfg["D_pass"],
                )
        result[name] = rec

    p = OUT / "a18_csv_measurements.json"
    p.write_text(json.dumps(result, ensure_ascii=False, indent=1), encoding="utf-8")

    print("\n=== CSV-обмеры владельца ===")
    print(f"{'образец':<18}{'вел.':<4}{'N':>4}{'среднее':>10}{'СКО':>9}"
          f"{'медиана':>10}{'MADσ':>9}{'мин':>9}{'макс':>9}")
    for name, rec in result.items():
        for kind in ("D", "P"):
            if not rec.get(kind):
                continue
            for tag in ("all", "clean"):
                s = rec[kind][tag]
                mark = "" if tag == "all" else " (роб.)"
                print(f"{name:<18}{kind + mark:<4}{s['n']:>4}{s['mean']:>10.3f}"
                      f"{s['std']:>9.3f}{s['median']:>10.3f}{s['mad_sigma']:>9.3f}"
                      f"{s['min']:>9.2f}{s['max']:>9.2f}")
    print("\n=== D против паспорта (гипотеза «+0.2 мкм») ===")
    for name, rec in result.items():
        if "D_vs_passport" not in rec:
            print(f"{name:<18} нет _D.csv")
            continue
        d = rec["D_vs_passport"]["clean"]
        print(f"{name:<18} паспорт {SAMPLES[name]['D_pass']:>5.1f} | "
              f"Δмед {d['delta_median']:+.3f} мкм | отн. {d['ratio_median']:.4f} | "
              f"Δср {d['delta_mean']:+.3f}")
    print(f"\nартефакт: {p.relative_to(REPO)}")
    return result


def _sample_frames(name):
    cfg = SAMPLES[name]
    raw = REPO / "data_pool" / name / "microscopy" / "raw"
    return sorted(raw.glob(cfg["raw_glob"]))


def cmd_sample(args):
    name = args.sample
    cfg = SAMPLES[name]
    frames = _sample_frames(name)
    if args.limit:
        frames = frames[:args.limit]
    print(f"\n=== {name}: {len(frames)} кадров, паспорт D={cfg['D_pass']} P={cfg['P_pass']} ===")
    per_frame, pooled_P, pooled_dev, pooled_W, fsum = [], [], [], [], []
    P_frame_means, P_frame_se = [], []
    for f in frames:
        r = process_frame(f, cfg["P_pass"], row_step=args.row_step, verbose=True)
        lv = r["levels"].get("0.50", {})
        ok_frame = (lv.get("n_wires", 0) >= 4 and lv.get("n_bad_wires", 99) == 0
                    and lv.get("n_narrow_inner", 99) == 0)
        if lv.get("n_wires", 0) >= 3:
            nper = len(lv["periods"])
            sdP = lv["period_stats"]["std"]
            fsum.append(dict(file=r["file"], n_wires=lv["n_wires"],
                             n_bad_wires=lv["n_bad_wires"], reliable=bool(ok_frame),
                             P_mean_um=lv["P_mean_um"],
                             # телескопическая СКО среднего периода: среднее шагов
                             # равно (x_N-x_1)/(N-1), ошибку дают только два конца
                             P_sem_um=sdP / max(nper, 1),
                             P_sem_naive_um=sdP / np.sqrt(nper),
                             P_std_um=sdP,
                             P_lsq_um=lv["P_lsq_um"],
                             dev_std_um=lv["dev_lsq"]["std"],
                             dev_owner_std_um=lv["dev_owner"]["std"],
                             dev_max_abs_um=lv["dev_owner"]["max_abs"],
                             ratio_std_over_P=lv["dev_lsq"]["ratio_std_over_P"],
                             acf1_periods=lv["acf1_periods"],
                             sigma_P_over_sigma_dev=lv["sigma_P_over_sigma_dev"],
                             trend_slope_um=lv["dev_trend"]["slope_um_per_wire"],
                             trend_t=lv["dev_trend"]["t"],
                             D_median_um=lv["width_stats"]["median"],
                             P_centroid_um=lv["centroid_control"]["P_mean_um"],
                             P_density_um=lv["P_density_um"],
                             n_units=lv["n_units"], n_merged=lv["n_merged"],
                             tilt_deg=r["slope"]["tilt_deg"],
                             wire=r["polarity"]["wire"]))
            if ok_frame:
                pooled_P += lv["periods"]
                pooled_dev += lv["dev_lsq_um"]
                pooled_W += lv["widths_um"]
                P_frame_means.append(lv["P_mean_um"])
                P_frame_se.append(sdP / max(nper, 1))
        r.pop("_prof"), r.pop("_pol"), r.pop("_p_approx")
        per_frame.append(r)

    # порог границы: как двигаются P и D (проверка оговорки ORCH)
    thr_scan = {}
    for key in [f"{lv:.2f}" for lv in LEVELS]:
        Ps, Ds, Ns = [], [], []
        for r in per_frame:
            v = r["levels"].get(key, {})
            if (v.get("n_wires", 0) >= 4 and v.get("n_bad_wires", 99) == 0
                    and v.get("n_narrow_inner", 99) == 0):
                Ps += v["periods"]
                Ds += v["widths_um"]
                Ns.append(v["n_wires"])
        if Ps:
            thr_scan[key] = dict(P=robust_stats(Ps, "P"), D=robust_stats(Ds, "D"),
                                 n_wires_total=int(sum(Ns)))
    keys = sorted(k for k in thr_scan if thr_scan[k])
    thr_effect = None
    if len(keys) >= 2:
        Pv = [thr_scan[k]["P"]["mean"] for k in keys]
        Dv = [thr_scan[k]["D"]["median"] for k in keys]
        thr_effect = dict(levels=[float(k) for k in keys],
                          P_mean=Pv, D_median=Dv,
                          dP_um=float(max(Pv) - min(Pv)),
                          dD_um=float(max(Dv) - min(Dv)),
                          dP_rel=float((max(Pv) - min(Pv)) / np.mean(Pv)),
                          dD_rel=float((max(Dv) - min(Dv)) / np.mean(Dv)))
        if "0.20" in thr_scan and "0.80" in thr_scan and thr_scan["0.20"] and thr_scan["0.80"]:
            # ширина фронта 20->80 % на ОДНУ границу: прямая мера размытия края
            thr_effect["edge_rise_2080_um"] = float(
                0.5 * (thr_scan["0.20"]["D"]["median"] - thr_scan["0.80"]["D"]["median"]))

    pooled_P = np.asarray(pooled_P)
    pooled_dev = np.asarray(pooled_dev)
    pooled_W = np.asarray(pooled_W)
    nmed = int(np.median([f["n_wires"] for f in fsum if f["reliable"]])) if any(
        f["reliable"] for f in fsum) else 0

    # средний период по образцу: взвешенное среднее кадровых средних
    P_best = P_best_se = None
    if P_frame_means:
        w = 1.0 / np.maximum(np.asarray(P_frame_se), 1e-9) ** 2
        P_best = float(np.sum(w * np.asarray(P_frame_means)) / np.sum(w))
        P_best_se = float(np.sqrt(1.0 / np.sum(w)))
        if len(P_frame_means) > 1:
            # разброс кадровых средних как независимая проверка внутренней ошибки
            P_best_se = max(P_best_se,
                            float(np.std(P_frame_means, ddof=1) / np.sqrt(len(P_frame_means))))

    # --- перекрывающиеся кадры -------------------------------------------------------
    # Кадры одного образца могут быть СДВИНУТЫМИ СНИМКАМИ ОДНОГО поля (микроскоп
    # подвинули на пол-кадра). Такие кадры не независимы: их нельзя считать за два
    # при оценке погрешности среднего. Ищем совпадающие цепочки периодов.
    seqs = [(r["file"], np.asarray(r["levels"]["0.50"]["periods"]))
            for r in per_frame if r["levels"].get("0.50", {}).get("n_wires", 0) >= 4]
    overlaps, dup_idx = [], set()
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i][1], seqs[j][1]
            best = 0
            for sh in range(-(len(a) - 1), len(b)):
                run = 0
                for k in range(max(0, -sh), min(len(a), len(b) - sh)):
                    if abs(a[k] - b[k + sh]) < 0.3:
                        run += 1
                    else:
                        best = max(best, run)
                        run = 0
                best = max(best, run)
            if best >= 4:
                overlaps.append(dict(a=seqs[i][0], b=seqs[j][0], matched=int(best)))
                dup_idx.add(j)
    n_independent = len(seqs) - len(dup_idx)

    # средний период по плотности проволок — включает и кадры со слипшимися
    # проволоками (их кратность учтена), поэтому не смещён отбраковкой коротких периодов
    dens = [f["P_density_um"] for f in fsum
            if f["P_density_um"] and f["n_units"] >= 4]
    P_density = dict(n_frames=len(dens),
                     mean=float(np.mean(dens)) if dens else None,
                     sem=(float(np.std(dens, ddof=1) / np.sqrt(len(dens)))
                          if len(dens) > 1 else None),
                     median=float(np.median(dens)) if dens else None,
                     values=dens)

    acf_pool = None
    if any(f["reliable"] for f in fsum):
        a = [f["acf1_periods"] for f in fsum if f["reliable"] and np.isfinite(f["acf1_periods"])]
        s = [f["sigma_P_over_sigma_dev"] for f in fsum
             if f["reliable"] and np.isfinite(f["sigma_P_over_sigma_dev"])]
        acf_pool = dict(acf1_periods_mean=float(np.mean(a)) if a else None,
                        acf1_periods_values=a,
                        sigma_ratio_mean=float(np.mean(s)) if s else None,
                        sigma_ratio_values=s)

    agg = dict(
        sample=name, calibration_um_per_px=CAL_UM_PER_PX,
        passport=dict(D_um=cfg["D_pass"], P_um=cfg["P_pass"]),
        expect_wire=cfg["expect_wire"],
        n_frames=len(frames),
        n_frames_reliable=int(sum(1 for f in fsum if f["reliable"])),
        frames_summary=fsum,
        P_best_um=P_best, P_best_se_um=P_best_se,
        P_density=P_density,
        overlapping_frames=overlaps, n_frames_independent=int(n_independent),
        pooled=dict(periods_um=pooled_P.tolist(), dev_lsq_um=pooled_dev.tolist(),
                    widths_um=pooled_W.tolist(),
                    P=robust_stats(pooled_P, "P_um"),
                    D=robust_stats(pooled_W, "D_um"),
                    dev=robust_stats(pooled_dev, "dev_um"),
                    sigma_P_over_P=(float(np.std(pooled_P, ddof=1) / np.mean(pooled_P))
                                    if pooled_P.size > 1 else None),
                    dev_std_um=float(np.std(pooled_dev, ddof=1)) if pooled_dev.size > 1 else None,
                    dev_max_abs_um=float(np.abs(pooled_dev).max()) if pooled_dev.size else None,
                    ratio_std_over_P=(float(np.std(pooled_dev, ddof=1) / np.mean(pooled_P))
                                      if pooled_dev.size > 1 else None),
                    duty_D_over_P=(float(np.median(pooled_W) / np.median(pooled_P))
                                   if pooled_W.size else None)),
        threshold_scan=thr_scan,
        threshold_effect=thr_effect,
        structure=acf_pool,
        null_model_iid=null_model_iid(max(nmed, 4)) if nmed else None,
        frames=per_frame,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"a18_{name}.json").write_text(
        json.dumps(agg, ensure_ascii=False, indent=1), encoding="utf-8")
    if fsum:
        plot_sample(agg, OUT / f"a18_{name}.png", title=f"A18 — {name}")
    # диагностический рисунок по первому кадру
    if frames:
        r0 = process_frame(frames[0], cfg["P_pass"], row_step=args.row_step)
        plot_frame(r0, OUT / f"a18_{name}_frame0.png", title=f"{name} / {frames[0].name}")

    print(f"\n  надёжных кадров {agg['n_frames_reliable']}/{len(frames)}, "
          f"независимых полей {n_independent}")
    for o in overlaps:
        print(f"    ⚠ перекрытие: {Path(o['a']).name} ~ {Path(o['b']).name} "
              f"({o['matched']} совпавших периодов подряд)")
    if P_best is not None:
        print(f"  P (взвеш. по чистым кадрам) = {P_best:.3f} ± {P_best_se:.3f} мкм")
    if P_density["mean"] is not None:
        print(f"  P (по плотности, {P_density['n_frames']} кадров) = {P_density['mean']:.3f}"
              + (f" ± {P_density['sem']:.3f}" if P_density["sem"] else "") + " мкм")
    if thr_effect:
        print(f"  порог f=0.2..0.8: ΔP = {thr_effect['dP_um']:.4f} мкм "
              f"({100*thr_effect['dP_rel']:.3f} %), ΔD = {thr_effect['dD_um']:.3f} мкм "
              f"({100*thr_effect['dD_rel']:.1f} %)")
    if pooled_P.size == 0:
        print("  НЕТ надёжных кадров — пул пуст")
        return agg
    print(f"  пул: N_периодов={pooled_P.size}  P_мед={np.median(pooled_P):.3f} мкм  "
          f"MADσ={1.4826*np.median(np.abs(pooled_P-np.median(pooled_P))):.3f}  "
          f"σ_P/P={100*np.std(pooled_P, ddof=1)/np.mean(pooled_P):.1f} %")
    if pooled_dev.size > 1:
        print(f"       СКО откл.={np.std(pooled_dev, ddof=1)*1000:.0f} нм  "
              f"СКО/P={100*np.std(pooled_dev, ddof=1)/np.mean(pooled_P):.2f} %")
    if pooled_W.size:
        print(f"       D_мед={np.median(pooled_W):.3f} мкм  "
              f"D/P={np.median(pooled_W)/np.median(pooled_P):.3f}")
    print(f"  артефакт: research/results/a18_microscopy/a18_{name}.json")
    return agg


def cmd_all(args):
    cmd_csv(args)
    for name in SAMPLES:
        a = argparse.Namespace(sample=name, limit=args.limit, row_step=args.row_step)
        cmd_sample(a)


def cmd_frame(args):
    p = Path(args.path)
    if not p.is_absolute():
        p = REPO / p
    name = args.sample
    cfg = SAMPLES[name]
    r = process_frame(p, cfg["P_pass"], row_step=args.row_step, verbose=True)
    for k, v in sorted(r["levels"].items()):
        if v.get("n_wires", 0) >= 3:
            print(f"  f={k}  N={v['n_wires']:2d}  P={v['P_mean_um']:.4f}  "
                  f"D_med={v['width_stats']['median']:.4f}  "
                  f"СКО_откл(МНК)={v['dev_lsq']['std']*1000:.0f} нм")
    plot_frame(r, OUT / f"a18_frame_{p.stem}.png", title=p.name)


def cmd_focus(args):
    """Фокус-серии test_grid_40_20 — зависимость D и P от положения фокуса."""
    name = "test_grid_40_20"
    cfg = SAMPLES[name]
    root = REPO / "data_pool" / name / "microscopy" / "raw" / "focus_series"
    series = sorted([d for d in root.iterdir() if d.is_dir()])
    res = {}
    for s in series:
        frames = sorted(s.glob("*.bmp"))
        rows = []
        print(f"\n--- {s.name}: {len(frames)} кадров")
        for f in frames:
            r = process_frame(f, cfg["P_pass"], row_step=args.row_step,
                              levels=(0.2, 0.5, 0.8), verbose=True)
            lv = r["levels"].get("0.50", {})
            l2, l8 = r["levels"].get("0.20", {}), r["levels"].get("0.80", {})
            if lv.get("n_wires", 0) >= 3:
                rise = None
                if l2.get("n_wires") == l8.get("n_wires") == lv["n_wires"]:
                    rise = 0.5 * (l2["width_stats"]["median"] - l8["width_stats"]["median"])
                rows.append(dict(file=f.name, n_wires=lv["n_wires"],
                                 n_bad_wires=lv["n_bad_wires"],
                                 n_units=lv["n_units"], n_merged=lv["n_merged"],
                                 P_mean_um=lv["P_mean_um"],
                                 P_density_um=lv["P_density_um"],
                                 D_median_um=lv["width_stats"]["median"],
                                 D20_median_um=(l2["width_stats"]["median"]
                                                if l2.get("n_wires", 0) >= 3 else None),
                                 D80_median_um=(l8["width_stats"]["median"]
                                                if l8.get("n_wires", 0) >= 3 else None),
                                 edge_rise_2080_um=rise,
                                 dev_std_um=lv["dev_lsq"]["std"],
                                 wire=r["polarity"]["wire"],
                                 contrast=r["slope"]["contrast"],
                                 tilt_deg=r["slope"]["tilt_deg"]))
        res[s.name] = rows
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "a18_test_grid_40_20_focus.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")

    plt = _mpl()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    allrise, allD, allD20, allP = [], [], [], []
    for sname, rows in res.items():
        rr = [r for r in rows if r["edge_rise_2080_um"] is not None]
        if not rr:
            continue
        x = [r["edge_rise_2080_um"] for r in rr]
        ax[0].plot(x, [r["D_median_um"] for r in rr], "o", label=sname, ms=5)
        ax[1].plot(x, [r["D20_median_um"] for r in rr], "s", label=sname, ms=5)
        ax[2].plot(x, [r["P_density_um"] for r in rr], "^", label=sname, ms=5)
        allrise += x
        allD += [r["D_median_um"] for r in rr]
        allD20 += [r["D20_median_um"] for r in rr]
        allP += [r["P_density_um"] for r in rr]
    for a, t, yl in zip(ax, ["D при f=0.50 против размытия края",
                             "D при f=0.20 (конвенция владельца)",
                             "P (по плотности) против размытия"],
                        ["D, мкм", "D, мкм", "P, мкм"]):
        a.set_xlabel("ширина фронта 20-80 %, мкм (размытие)")
        a.set_ylabel(yl)
        a.set_title(t, fontsize=10)
    ax[0].legend(fontsize=7)
    for a, y in zip(ax[:2], [20.0, 20.0]):
        a.axhline(y, color="k", ls=":", lw=1)
    fig.tight_layout()
    fig.savefig(OUT / "a18_test_grid_40_20_focus.png", dpi=110)
    plt.close(fig)

    def _corr(x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 3:
            return None
        return dict(n=int(m.sum()), r=float(np.corrcoef(x[m], y[m])[0, 1]),
                    slope=float(np.polyfit(x[m], y[m], 1)[0]),
                    intercept=float(np.polyfit(x[m], y[m], 1)[1]))

    summary = dict(D50_vs_rise=_corr(allrise, allD),
                   D20_vs_rise=_corr(allrise, allD20),
                   P_vs_rise=_corr(allrise, allP),
                   rise=robust_stats(allrise, "rise"),
                   D50=robust_stats(allD, "D50"),
                   D20=robust_stats(allD20, "D20"),
                   P_density=robust_stats(allP, "P_density"))
    res["_summary"] = summary
    (OUT / "a18_test_grid_40_20_focus.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print("\n--- зависимость от размытия края (все фокус-серии) ---")
    for k in ("D50_vs_rise", "D20_vs_rise", "P_vs_rise"):
        c = summary[k]
        if c:
            print(f"  {k}: n={c['n']} r={c['r']:+.3f} наклон={c['slope']:+.3f} "
                  f"экстраполяция к нулевому размытию = {c['intercept']:.3f} мкм")
    print("\nартефакт: research/results/a18_microscopy/a18_test_grid_40_20_focus.json")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("csv");    p.set_defaults(fn=cmd_csv)
    p = sub.add_parser("sample"); p.add_argument("sample", choices=list(SAMPLES))
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--row-step", type=int, default=1)
    p.set_defaults(fn=cmd_sample)
    p = sub.add_parser("all")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--row-step", type=int, default=1)
    p.set_defaults(fn=cmd_all)
    p = sub.add_parser("frame")
    p.add_argument("path")
    p.add_argument("--sample", choices=list(SAMPLES), required=True)
    p.add_argument("--row-step", type=int, default=1)
    p.set_defaults(fn=cmd_frame)
    p = sub.add_parser("focus")
    p.add_argument("--row-step", type=int, default=2)
    p.set_defaults(fn=cmd_focus)

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
