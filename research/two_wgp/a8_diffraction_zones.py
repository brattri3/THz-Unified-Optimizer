"""A8 — оптическая дифракция от PureWave: сравнение восстановленной и заводской зон.

КОНТЕКСТ (метаданные владельца 2026-08-04). Поляризатор PureWave был ВОССТАНОВЛЕН:
центральная зона потеряла натяжение проводов (отслоение клея от оправы), провода
провисли. Владелец переклеил и деформацией оправы подтянул центральные провода;
периодичность восстановлена не полностью, по краям сохранились заводские участки
с хорошей периодичностью. Разница видна невооружённым глазом по дифракции лазерной
указки: снимки `good N.jpg` (заводская зона) и `bad N.jpg` (восстановленная).

ЗАЧЕМ ЭТО ЗДЕСЬ. Проволочная сетка — амплитудная дифракционная решётка. В видимом
свете (lambda = 0.65 мкм << P = 25.5 мкм) она даёт много порядков, и КАЧЕСТВО картины
прямо измеряет то, что в ТГц-модели Бланко заложено как ИДЕАЛЬНАЯ периодичность:
  * положение порядков   -> период P (и его локальные вариации);
  * резкость порядков    -> когерентная длина решётки / разброс периода;
  * контраст между порядками и диффузный фон -> доля непериодической структуры.
Для загадки D_eff (все шесть образцов требуют D_eff/D_phys = 0.37..0.59) это первый
случай, когда нерегулярность известна независимо и её можно измерить.

ЧТО СЧИТАЕТСЯ
-------------
1. Ось дифракционной линии — взвешенная по яркости регрессия (линия на снимках наклонена).
2. Профиль I(x) вдоль оси: суммирование поперёк полосы, вычитание локального фона.
3. Позиции порядков (`find_peaks` по prominence в логарифме) + шаг между ними.
4. Метрики качества, сравнимые между снимками:
   - `n_orders` — сколько порядков различимо;
   - `step_px`, `step_rel_std` — средний шаг и его относительный разброс (регулярность);
   - `acf_step_px` — период по автокорреляции (устойчив к пропущенным пикам);
   - `visibility` — средняя видность (I_max-I_min)/(I_max+I_min) между соседними порядками;
   - `diffuse_fraction` — доля энергии профиля ВНЕ пиков (гало от непериодичности);
   - `order_fwhm_px` — средняя ширина порядка (уширение = разброс периода).

ГЕОМЕТРИЯ (числа владельца): lambda = 640..660 нм (указка), L = 315 мм до экрана,
расстояние между двумя первыми порядками = 1.75 мм. Пересчёт d = m*lambda*L/x
делается в обе стороны и сверяется с микрофото (P = 25.5 +- 0.9 мкм) — см. вывод
`geometry_check`: он показывает, что 1.75 мм НЕ согласуется с микрофото, а 17.5 мм —
согласуется, и это надо перепроверить на установке.

ВАЖНО ПРО МАСШТАБ. В кадре нет линейки, поэтому абсолютный шаг в мм из снимка НЕ
восстановим — пиксельный шаг зависит от расстояния съёмки и зума. Поэтому все
абсолютные выводы о периоде идут от числа владельца, а снимки дают ОТНОСИТЕЛЬНЫЕ
величины: сравнение зон между собой и регулярность внутри зоны. Отношение шагов
good/bad — величина безразмерная и от масштаба не зависит только при одинаковых
условиях съёмки, поэтому оно приводится с явной оговоркой.

Зона A, только чтение снимков. Артефакт -> research/results/two_wgp/a8_diffraction.json

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a8_diffraction_zones.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.signal import find_peaks

REPO = Path(__file__).resolve().parents[2]
IMG_DIR = REPO / "research" / "validation" / "micrographs" / "purewave"
OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_NM = (640.0, 660.0)          # надпись на указке
L_MM = 315.0                        # поляризатор -> экран
# Расстояние между порядками m=+1 и m=-1. Владелец уточнил 2026-08-04: 17.5 мм
# (в первом сообщении было «1,75 мм» — описка, подтверждено им же после сверки
# с микрофото и с оценкой расстояния съёмки по кадру).
DX_FIRST_ORDERS_MM = 17.5
DX_ERR_MM = 0.5                     # реалистичная точность отсчёта по экрану
L_ERR_MM = 2.0
P_MICROGRAPH_UM = 25.5              # research/validation/MANIFEST.md
P_MICROGRAPH_ERR_UM = 0.9


# --------------------------------------------------------------- профиль
def line_profile(path, half_width=18, bg_offset=90):
    """Профиль интенсивности вдоль дифракционной линии.

    Возвращает (x, I, meta). Ось находится взвешенной регрессией по ярким пикселям,
    поэтому наклон снимка учитывается. Фон берётся из двух полос, отстоящих от оси
    на `bg_offset` пикселей, и вычитается — так убирается вуаль и засветка кадра.
    """
    im = Image.open(path)
    a = np.asarray(im).astype(np.float64)
    # лазер красный: красный канал даёт лучший контраст; при насыщении центра это
    # не мешает, положения БОКОВЫХ порядков определяются не по нему.
    I = a[:, :, 0] if a.ndim == 3 else a
    h, w = I.shape

    thr = np.percentile(I, 99.9)
    ys, xs = np.nonzero(I >= thr)
    wts = I[ys, xs]
    # взвешенная прямая y = k*x + b
    k, b = np.polyfit(xs, ys, 1, w=wts)

    xg = np.arange(w)
    yc = k * xg + b
    rows = np.arange(h)[:, None]
    dy = rows - yc[None, :]
    core = np.abs(dy) <= half_width
    bg_band = (np.abs(dy) >= bg_offset) & (np.abs(dy) <= bg_offset + 40)

    prof = np.where(core, I, 0.0).sum(axis=0)
    n_core = core.sum(axis=0)
    bg = np.where(bg_band, I, np.nan)
    with np.errstate(invalid="ignore"):
        bg_med = np.nanmedian(bg, axis=0)
    bg_med = np.where(np.isfinite(bg_med), bg_med, 0.0)
    prof = prof - bg_med * n_core

    meta = {"file": Path(path).name, "shape": [int(h), int(w)],
            "line_slope": float(k), "line_intercept_px": float(b),
            "line_tilt_deg": float(np.degrees(np.arctan(k))),
            "saturated_px": int(np.sum(I >= 254.0)),
            "half_width_px": half_width}
    return xg.astype(float), prof, meta


# ДВА МАСШТАБА В КАДРЕ. Визуальный разбор выпрямленной полосы (см. артефакт) показывает
# структуру двух масштабов: мелкую ~20 px и крупную ~100 px. Мелкая присутствует и в
# заводской, и в восстановленной зоне одинаково и не образует эквидистантного ряда —
# это СПЕКЛ (лазерная указка многомодовая, экран диффузный). Крупная — разделённые
# тёмными промежутками «острова», эквидистантные у заводской зоны, — это и есть
# дифракционные порядки. Поэтому анализ ведётся на масштабе 40..250 px, а спекл
# подавляется сглаживанием.
SMOOTH_PX = 15          # подавление спекла
ENVELOPE_PX = 301       # огибающая (спад яркости от центра) — её делим
PERIOD_MIN_PX = 40
PERIOD_MAX_PX = 250


def modulation(prof, i0, half_window=600):
    """Модуляция профиля = сглаженный профиль / медленная огибающая.

    Деление на огибающую убирает спад яркости от центра к краям (он не несёт информации
    о регулярности) и делает `good` и `bad` сравнимыми при разной экспозиции.
    Сглаживание на SMOOTH_PX убирает спекл, оставляя структуру порядков.
    """
    from scipy.ndimage import median_filter, uniform_filter1d
    lo, hi = max(i0 - half_window, 0), min(i0 + half_window, len(prof))
    seg = np.asarray(prof[lo:hi], float)
    sm = uniform_filter1d(seg, size=SMOOTH_PX, mode="nearest")
    env = median_filter(sm, size=ENVELOPE_PX, mode="nearest")
    env = np.maximum(env, max(1.0, 0.01 * np.median(np.abs(sm))))
    return sm / env, lo, hi


def spectral_regularity(mod):
    """Насколько модуляция ПЕРИОДИЧНА: спектр мощности и доля энергии в главном пике.

    Идеальная решётка -> узкий пик на частоте порядков (мощность сосредоточена).
    Разрушенная периодичность -> энергия размазана (гало от непериодической структуры).
    Возвращает период главного пика (px), долю энергии в пике +-2 бина, отношение
    пика к медиане спектра и глубину модуляции.
    """
    m = mod - mod.mean()
    n = len(m)
    S = np.abs(np.fft.rfft(m * np.hanning(n))) ** 2
    freq = np.fft.rfftfreq(n, d=1.0)
    ok = (freq > 1.0 / PERIOD_MAX_PX) & (freq < 1.0 / PERIOD_MIN_PX)
    if not np.any(ok):
        return None
    idx = np.where(ok)[0]
    j = idx[int(np.argmax(S[idx]))]
    lo, hi = max(j - 2, idx[0]), min(j + 3, idx[-1] + 1)   # не выходим за рабочий диапазон
    return {
        "period_px": float(1.0 / freq[j]),
        "peak_over_median": float(S[j] / np.median(S[idx])),
        "energy_fraction_in_peak": float(S[lo:hi].sum() / S[idx].sum()),
        "modulation_rms": float(np.std(mod)),
    }


def analyse_image(path, verbose=True):
    x, prof, meta = line_profile(path)
    if prof.max() - np.median(prof) <= 0:
        return {**meta, "error": "пустой профиль"}

    i0 = int(np.argmax(prof))
    mod, lo, hi = modulation(prof, i0)
    spec = spectral_regularity(mod)

    # порядки ищем в МОДУЛЯЦИИ (в исходном профиле их скрывает спад яркости);
    # distance = половина ожидаемого периода, чтобы не считать спекл за порядок
    sm = mod
    peaks, props = find_peaks(sm, prominence=0.05, distance=PERIOD_MIN_PX // 2)

    d = np.diff(peaks).astype(float)
    step = float(np.median(d)) if len(d) else None
    step_std = float(np.std(d, ddof=1)) if len(d) > 1 else None

    # видность между соседними порядками (по модуляции — не зависит от экспозиции)
    vis = []
    for a_, b_ in zip(peaks[:-1], peaks[1:]):
        if b_ - a_ < 3:
            continue
        imin, imax = sm[a_:b_].min(), min(sm[a_], sm[b_])
        if imax + imin > 0:
            vis.append((imax - imin) / (imax + imin))
    visibility = float(np.mean(vis)) if vis else None

    # FWHM боковых порядков в модуляции (нулевой исключён: насыщен)
    fwhms = []
    for pk in peaks:
        if abs(pk - (i0 - lo)) < 3:
            continue
        half = (sm[pk] + 1.0) / 2.0
        l_ = pk
        while l_ > 0 and sm[l_] > half:
            l_ -= 1
        r_ = pk
        while r_ < len(sm) - 1 and sm[r_] > half:
            r_ += 1
        fwhms.append(r_ - l_)
    fwhm = float(np.median(fwhms)) if fwhms else None

    lgmax = float(np.log10(max(prof.max(), 1.0)))
    res = {
        **meta,
        "zone": "заводская (good)" if Path(path).stem.startswith("good") else "восстановленная (bad)",
        "n_orders": int(len(peaks)),
        "order0_index_px": i0,
        "window_px": [int(lo), int(hi)],
        "step_px": step,
        "step_std_px": step_std,
        "step_rel_std": (step_std / step) if (step and step_std) else None,
        "spectral": spec,
        "visibility": visibility,
        "order_fwhm_px": fwhm,
        "fwhm_over_step": (fwhm / step) if (fwhm and step) else None,
        "peak_positions_px": [int(v + lo) for v in peaks],
        "dynamic_range_log10": lgmax - float(np.log10(max(np.median(prof), 1.0))),
    }
    if verbose:
        s = spec or {}
        print(f"  {res['file']:<12} [{res['zone']:<22}] порядков {res['n_orders']:>3}, "
              f"шаг {step if step else float('nan'):>6.2f} px "
              f"(разброс {100 * (res['step_rel_std'] or 0):>4.1f}%), "
              f"период по спектру {s.get('period_px', float('nan')):>6.2f} px, "
              f"пик/медиана {s.get('peak_over_median', float('nan')):>7.1f}, "
              f"глубина модуляции {s.get('modulation_rms', float('nan')):.3f}, "
              f"видность {visibility if visibility else float('nan'):.3f}")
    return res


def reproducibility(files, prefix, max_shift=60):
    """Попарная корреляция профилей внутри зоны — воспроизводима ли картина.

    Устойчивая решётка даёт ОДНУ И ТУ ЖЕ картину при повторной съёмке того же участка;
    разрушенная — картину, которая меняется от кадра к кадру (гало формируется случайной
    непериодической структурой). Сдвиг подбирается, масштаб — нет, поэтому число
    занижено для кадров с разным зумом; это указано в артефакте.
    """
    mods = []
    for p in files:
        if not p.stem.startswith(prefix):
            continue
        _x, prof, _m = line_profile(p)
        i0 = int(np.argmax(prof))
        mod, _lo, _hi = modulation(prof, i0)
        mods.append((p.stem, mod - mod.mean()))
    pairs = []
    for i in range(len(mods)):
        for j in range(i + 1, len(mods)):
            a, b = mods[i][1], mods[j][1]
            n = min(len(a), len(b))
            a, b = a[:n], b[:n]
            best = -1.0
            for s in range(-max_shift, max_shift + 1):
                aa = a[max(s, 0): n + min(s, 0)]
                bb = b[max(-s, 0): n + min(-s, 0)]
                if len(aa) < 100:
                    continue
                den = np.linalg.norm(aa) * np.linalg.norm(bb)
                if den > 0:
                    best = max(best, float(aa @ bb / den))
            pairs.append({"a": mods[i][0], "b": mods[j][0], "max_corr": best})
    return {"pairs": pairs,
            "mean_corr": float(np.mean([p["max_corr"] for p in pairs])) if pairs else None}


def make_figure(files, out_png):
    """Профили модуляции: заводская зона vs восстановленная, на одной шкале."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, sharey=True)
    for ax, prefix, title in ((axes[0], "good", "Заводская зона (good)"),
                              (axes[1], "bad", "Восстановленная зона (bad)")):
        for p in files:
            if not p.stem.startswith(prefix):
                continue
            x, prof, _m = line_profile(p)
            i0 = int(np.argmax(prof))
            mod, lo, hi = modulation(prof, i0)
            ax.plot(np.arange(lo, hi) - i0, mod, lw=0.9, label=p.stem)
        ax.set_title(title, fontsize=10)
        ax.axhline(1.0, color="k", lw=0.5, alpha=0.4)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=3)
        ax.set_ylabel("I / огибающая")
    axes[1].set_xlabel("положение вдоль дифракционной линии, px (0 = нулевой порядок)")
    fig.suptitle("PureWave: дифракция лазерной указки, модуляция профиля\n"
                 "(глубокая регулярная модуляция = сохранившаяся периодичность решётки)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# ------------------------------------------------------------- геометрия
def geometry_check(step_px_measured=None):
    """Пересчёт периода из числа владельца и сверка с микрофото.

    Малые углы: sin(theta) ~ x/L, условие максимума d*sin(theta_m) = m*lambda.
      * если 1.75 мм — расстояние между СОСЕДНИМИ порядками: d = lambda*L/x;
      * если между m=+1 и m=-1 (т.е. полушаг = x/2):        d = 2*lambda*L/x.
    """
    lam = np.mean(LAMBDA_NM) * 1e-6                       # мм
    lam_err = (LAMBDA_NM[1] - LAMBDA_NM[0]) / 2 * 1e-6    # надпись 640..660 -> +-10 нм
    out = {"lambda_mm": lam, "L_mm": L_MM,
           "micrograph_P_um": P_MICROGRAPH_UM, "micrograph_P_err_um": P_MICROGRAPH_ERR_UM,
           "variants": []}

    # Основной результат: 17.5 мм между m=+1 и m=-1 => полушаг x/2, d = 2*lambda*L/x.
    d_um = 2 * lam * L_MM / DX_FIRST_ORDERS_MM * 1000.0
    rel = np.sqrt((lam_err / lam) ** 2 + (L_ERR_MM / L_MM) ** 2
                  + (DX_ERR_MM / DX_FIRST_ORDERS_MM) ** 2)
    d_err = d_um * rel
    diff = P_MICROGRAPH_UM - d_um
    sig = float(np.hypot(d_err, P_MICROGRAPH_ERR_UM))
    out["primary"] = {
        "interpretation": "17.5 мм между m=+1 и m=-1 (уточнено владельцем 2026-08-04)",
        "P_diffraction_um": d_um, "P_diffraction_err_um": d_err,
        "rel_err": rel,
        "error_budget": {"lambda": lam_err / lam, "L": L_ERR_MM / L_MM,
                         "x": DX_ERR_MM / DX_FIRST_ORDERS_MM},
        "micrograph_minus_diffraction_um": diff,
        "discrepancy_sigma": diff / sig,
        "consistent": bool(abs(diff) <= 2 * sig),
        "max_order": int(np.floor(d_um * 1e-3 / lam)),
        "caveat": "оптический период относится К ТОЙ ЗОНЕ, куда светила указка; микрофото — "
                  "к другому (неизвестному) участку. При пространственной неоднородности "
                  "образца это РАЗНЫЕ величины, а не два измерения одной.",
    }

    for label, x_mm, factor in (
            ("уточнено: 17.5 мм между m=+1 и m=-1", DX_FIRST_ORDERS_MM, 2.0),
            ("если бы это был шаг СОСЕДНИХ порядков", DX_FIRST_ORDERS_MM, 1.0),
            ("первоначальная описка: 1.75 мм между m=+1 и m=-1", DX_FIRST_ORDERS_MM / 10, 2.0),
            ("первоначальная описка: 1.75 мм между СОСЕДНИМИ", DX_FIRST_ORDERS_MM / 10, 1.0)):
        d_mm = factor * lam * L_MM / x_mm
        d_um = d_mm * 1000.0
        out["variants"].append({
            "interpretation": label, "x_mm": x_mm,
            "d_um": d_um, "ratio_to_micrograph": d_um / P_MICROGRAPH_UM,
            "consistent": bool(abs(d_um - P_MICROGRAPH_UM) <= 3 * P_MICROGRAPH_ERR_UM),
        })
    # Проверка правдоподобия по кадру: если шаг соседних порядков равен `step_px`
    # пикселей и одновременно `s_mm` миллиметрам, то весь кадр (4000 px) покрывает
    # 4000/(step_px/s_mm) мм. Для телефона горизонтальное поле зрения примерно
    # 1.4*расстояние (~70 град), отсюда оценка расстояния съёмки. Она грубая, но
    # различает варианты, отличающиеся в 10 раз.
    if step_px_measured:
        out["frame_plausibility"] = []
        out["frame_plausibility_caveat"] = (
            "оценка расстояния предполагает поле зрения ~70° (без зума); при съёмке с зумом "
            "поле уже, и то же расстояние даёт больший масштаб. Поэтому проверка различает "
            "варианты, отличающиеся в 10 раз, но не в 2")
        for label, s_mm in (("шаг 8.75 мм (уточнённые 17.5 между ±1)", DX_FIRST_ORDERS_MM / 2),
                            ("шаг 8.03 мм (из микрофото P=25.5)",
                             lam * L_MM / (P_MICROGRAPH_UM * 1e-3)),
                            ("шаг 0.875 мм (описка 1.75 между ±1)", DX_FIRST_ORDERS_MM / 20),
                            ("шаг 1.75 мм (описка, как соседние)", DX_FIRST_ORDERS_MM / 10)):
            px_per_mm = step_px_measured / s_mm
            frame_mm = 4000.0 / px_per_mm
            out["frame_plausibility"].append({
                "assumption": label, "step_mm": s_mm,
                "px_per_mm": px_per_mm, "frame_width_mm": frame_mm,
                "implied_camera_distance_mm": frame_mm / 1.4,
            })

    # что ОЖИДАЕТСЯ увидеть при P по микрофото
    out["expected_from_micrograph"] = {
        "step_adjacent_mm": lam * L_MM / (P_MICROGRAPH_UM * 1e-3),
        "span_pm1_mm": 2 * lam * L_MM / (P_MICROGRAPH_UM * 1e-3),
        "max_order": int(np.floor(P_MICROGRAPH_UM * 1e-3 / lam)),
        "comment": "при P=25.5 мкм и L=315 мм соседние порядки отстоят на ~8 мм, "
                   "т.е. +1 и -1 разнесены на ~16 мм",
    }
    return out


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=str(IMG_DIR))
    ap.add_argument("--out", default=str(OUT / "a8_diffraction.json"))
    args = ap.parse_args()

    d = Path(args.dir)
    files = sorted([p for p in d.glob("*.jpg")
                    if p.stem.startswith(("good", "bad"))])
    if not files:
        print(f"ОСТАНОВ: в {d} нет снимков good*/bad*.jpg")
        return 1

    print(f"Снимки дифракции: {len(files)} шт. в {d.relative_to(REPO)}\n")
    rows = [analyse_image(p) for p in files]

    good = [r for r in rows if r.get("zone", "").startswith("заводская") and "step_px" in r]
    bad = [r for r in rows if r.get("zone", "").startswith("восстановленная") and "step_px" in r]

    def agg(rs, key):
        v = [r[key] for r in rs if r.get(key) is not None]
        return {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                "n": len(v)} if v else None

    # спектральные метрики лежат во вложенном словаре — поднимаем на верхний уровень
    for r in rows:
        for k, v in (r.get("spectral") or {}).items():
            r[f"spec_{k}"] = v

    comparison = {}
    for key in ("n_orders", "step_px", "step_rel_std", "visibility",
                "fwhm_over_step", "dynamic_range_log10",
                "spec_period_px", "spec_peak_over_median",
                "spec_energy_fraction_in_peak", "spec_modulation_rms"):
        g, b = agg(good, key), agg(bad, key)
        comparison[key] = {"good": g, "bad": b,
                           "bad_over_good": (b["mean"] / g["mean"]
                                             if g and b and g["mean"] else None)}

    repro = {"good": reproducibility(files, "good"), "bad": reproducibility(files, "bad")}
    print(f"\n  воспроизводимость картины между снимками одной зоны "
          f"(макс. корреляция профилей): good {repro['good']['mean_corr']:+.3f}, "
          f"bad {repro['bad']['mean_corr']:+.3f}")

    steps = [r["spectral"]["period_px"] for r in rows
             if r.get("spectral") and r["spectral"].get("period_px")]
    geo = geometry_check(float(np.median(steps)) if steps else None)

    print("\n--- СРАВНЕНИЕ ЗОН (заводская vs восстановленная) ---")
    for key, c in comparison.items():
        if not (c["good"] and c["bad"]):
            continue
        print(f"  {key:<20} good {c['good']['mean']:>8.4f} +- {c['good']['std']:<8.4f} "
              f"bad {c['bad']['mean']:>8.4f} +- {c['bad']['std']:<8.4f}"
              + (f"   отношение bad/good = {c['bad_over_good']:.3f}"
                 if c["bad_over_good"] else ""))

    print("\n--- ГЕОМЕТРИЯ: период из дифракции vs микрофото ---")
    pr = geo["primary"]
    print(f"  ОПТИЧЕСКИЙ ПЕРИОД: P = {pr['P_diffraction_um']:.2f} ± {pr['P_diffraction_err_um']:.2f} мкм "
          f"({pr['interpretation']})")
    print(f"    бюджет ошибки: λ {100 * pr['error_budget']['lambda']:.1f}%, "
          f"L {100 * pr['error_budget']['L']:.1f}%, отсчёт x {100 * pr['error_budget']['x']:.1f}% "
          f"⇒ {100 * pr['rel_err']:.1f}%")
    print(f"    микрофото − дифракция = {pr['micrograph_minus_diffraction_um']:+.2f} мкм "
          f"({pr['discrepancy_sigma']:+.2f}σ) — "
          f"{'совместимо' if pr['consistent'] else 'РАСХОЖДЕНИЕ'}; порядков до m = {pr['max_order']}")
    e = geo["expected_from_micrograph"]
    print(f"  при P = {P_MICROGRAPH_UM} мкм (микрофото) и L = {L_MM} мм ожидается: "
          f"соседние порядки {e['step_adjacent_mm']:.2f} мм, "
          f"размах +1..-1 {e['span_pm1_mm']:.2f} мм, порядков до m = {e['max_order']}")
    for v in geo["variants"]:
        print(f"  {v['interpretation']:<46} -> d = {v['d_um']:>8.1f} мкм "
              f"({v['ratio_to_micrograph']:>5.2f} x микрофото)"
              f"{'   <-- согласуется' if v['consistent'] else ''}")
    if geo.get("frame_plausibility"):
        print(f"\n  Проверка правдоподобия по кадру (медианный шаг "
              f"{np.median(steps):.1f} px, кадр 4000 px):")
        for fp in geo["frame_plausibility"]:
            print(f"    {fp['assumption']:<34} -> кадр {fp['frame_width_mm']:>7.0f} мм, "
                  f"съёмка с ~{fp['implied_camera_distance_mm'] / 10:>5.0f} см")

    payload = {"task": "A8_diffraction_zones", "session": "A", "nick": "A3",
               "images": rows, "comparison": comparison, "reproducibility": repro,
               "geometry_check": geo,
               "method_note": "для количественного вывода нужна пересъёмка с ЛИНЕЙКОЙ в "
                              "кадре и фиксированным расстоянием: сейчас масштаб кадров "
                              "различается (период главного пика 150..240 px), поэтому "
                              "сравнимы только безразмерные метрики",
               "scale_caveat": "в кадре нет линейки: пиксельный масштаб зависит от "
                               "расстояния и зума, поэтому абсолютный период из снимка "
                               "не восстанавливается; сравнимы только безразмерные метрики"}
    Path(args.out).write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                              encoding="utf-8")
    png = Path(args.out).with_suffix(".png")
    make_figure(files, png)
    print(f"\nАртефакт: {Path(args.out).relative_to(REPO)}")
    print(f"График:   {png.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
