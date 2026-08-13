# -*- coding: utf-8 -*-
u"""Бит-в-бит сверка порта фита с оригиналом сессии A.

``track_viewer/core/fit_malus.py`` — копия аппарата из
``research/two_wgp/malus_linearization.py``, а не импорт: комплект для
Windows 7 каталога ``research/`` не содержит и содержать не может. Копия
означает риск разъехаться, поэтому сверка обязана быть машинной.

Проверяется совпадение на ОДНИХ И ТЕХ ЖЕ входах — то есть что порт не изменил
метод. Согласие с данными проверяет приёмка (``track_viewer.selftest``), это
разные вопросы.

Скрипт работает только в репозитории. Внутри комплекта ``research/`` нет, и
запускать его там незачем — поэтому он лежит в ``tools/``, а не в приёмке.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe track_viewer\\tools\\verify_fit_malus.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "research", "two_wgp"))

from track_viewer.core import fit_malus as port          # noqa: E402
from track_viewer.core.compat import setup_console       # noqa: E402


def _synth(theta_deg, t_perp, t_par, theta0_deg):
    t = (np.asarray(theta_deg, float) - theta0_deg) * (np.pi / 180.0)
    E = t_perp * np.cos(t) ** 2 + t_par * np.sin(t) ** 2
    return np.abs(E) ** 2


#: Наборы углов: полный, урезанный производственный и заведомо неравномерный —
#: на последнем гармоники 2 и 4 неортогональны, и расхождение реализаций
#: вылезло бы именно там.
ANGLE_SETS = (
    ("полный -90..90 шаг 10", np.arange(-90.0, 91.0, 10.0)),
    ("урезанный 0..90", np.arange(0.0, 91.0, 10.0)),
    ("неравномерный", np.array([-100.0, -94.0, -90.0, -86.0, -80.0, -40.0,
                                0.0, 10.0, 55.0, 80.0, 88.0, 90.0, 92.0, 100.0])),
)

CASES = (
    (0.93 + 0.11j, 0.031 * np.exp(1j * 0.7), 1.35),
    (0.99 + 0.0j, 0.012 * np.exp(1j * 2.9), -0.463),
    (0.70 - 0.30j, 0.150 * np.exp(1j * 0.0), 12.7),
    (0.95 + 0.05j, 0.0005 * np.exp(1j * 1.1), -7.25),
)


def _max_diff(a, b):
    u"""Наибольшее расхождение по величинам, которые есть в ОБОИХ словарях.

    Порт сознательно у́же оригинала: часть производных величин в инструменте
    у прибора не нужна. Пропущенное собирается отдельно и печатается — молча
    сузить сверку значило бы сделать её бессмысленной.
    """
    worst = 0.0
    for k in set(a) & set(b):
        x, y = a[k], b[k]
        if isinstance(x, bool) or isinstance(y, bool):
            if x != y:
                return float("inf")
            continue
        if not (isinstance(x, float) and isinstance(y, float)):
            continue
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        worst = max(worst, abs(x - y))
    return worst


def main():
    setup_console()
    omitted = set()
    try:
        import malus_linearization as orig
    except ImportError as exc:
        print(u"оригинал сессии A недоступен (%s)" % exc)
        print(u"скрипт работает только в репозитории — в комплекте research/ нет")
        return 2

    print(u"сверка track_viewer/core/fit_malus.py с research/two_wgp/malus_linearization.py")
    print(u"")
    worst = 0.0
    checks = 0
    failed = 0

    for set_name, th in ANGLE_SETS:
        for t_perp, t_par, t0 in CASES:
            U = _synth(th, t_perp, t_par, t0)
            for weights in ("uniform", "relative"):
                sigma = None if weights == "uniform" else np.maximum(U, U.max() * 1e-12)
                for order in (2, 4):
                    if len(th) < (5 if order >= 4 else 3):
                        continue
                    p_a, cov_a, info_a = orig.fit_harmonics(th, U, sigma, order)
                    p_b, cov_b, info_b = port.fit_harmonics(th, U, sigma, order)
                    checks += 1
                    d = float(np.max(np.abs(p_a - p_b)))
                    dc = float(np.max(np.abs(cov_a - cov_b)))
                    if order >= 4:
                        ca = orig.channel_params(p_a, cov_a)
                        cb = port.channel_params(p_b, cov_b)
                        omitted.update(set(ca) - set(cb))
                        dp = _max_diff(ca, cb)
                        ta = orig.theta0_error(p_a, cov_a, info_a["chi2_red"])
                        tb = port.theta0_error(p_b, cov_b, info_b["chi2_red"])
                        dp = max(dp, abs(ta - tb) if np.isfinite(ta) else 0.0)
                    else:
                        ca = orig.malus_params(p_a, cov_a)
                        cb = port.malus_params(p_b, cov_b)
                        omitted.update(set(ca) - set(cb))
                        dp = _max_diff(ca, cb)
                    da = orig.design_diagnostics(th, order)
                    db = port.design_diagnostics(th, order)
                    dd = max(abs(da[k] - db[k]) for k in
                             ("cond_XtX", "sigma_theta0_deg", "max_overlap_h2_h4"))
                    m = max(d, dc, dp, dd)
                    worst = max(worst, m)
                    if m > 0.0:
                        failed += 1
                        print(u"[РАСХОЖДЕНИЕ] %s, order=%d, веса %s: %.3e"
                              % (set_name, order, weights, m))

    print(u"сверок: %d" % checks)
    print(u"максимальное расхождение по всем величинам: %.3e" % worst)
    if omitted:
        print(u"в порт не взяты (и потому не сверяются): %s"
              % u", ".join(sorted(omitted)))
    if failed:
        print(u"")
        print(u"ПОРТ РАЗОШЁЛСЯ С ОРИГИНАЛОМ — числа инструмента больше не совпадают")
        print(u"с числами сессии A, дальше идти нельзя")
        return 1
    print(u"")
    print(u"=== порт совпадает с оригиналом бит-в-бит ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
