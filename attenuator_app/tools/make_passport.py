"""Сборка паспорта прибора.

Запуск (из корня THz-Unified-Optimizer):
    .venv\\Scripts\\python.exe -m attenuator_app.tools.make_passport

ИСТОЧНИК ИСТИНЫ — официальный паспорт продукта
`THz-Spectroscopy-Python-Manual2/text/passport.tex` (ATT-11-16-CA85, S/N #02721).
Оттуда берутся значения параметров и их паспортные 95 % интервалы.

Файл `passports/_calib_raw.json` (наши фиты по шести сеансам) используется
ТОЛЬКО как источник консервативной оценки воспроизводимости: если разброс
между сеансами шире паспортного интервала, берётся он. Паспортный интервал
получен на отфильтрованном по дрейфу ансамбле и потому оптимистичен для
повседневной работы, где прибор переставляют и заново юстируют.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PASSPORTS = HERE.parent / "passports"
sys.path.insert(0, str(HERE.parents[1]))

from attenuator_app.core.passport import Passport, Param, ScaleSpec, film_for_aperture

# --- ПАСПОРТ ПРОДУКТА (passport.tex, таблица разд. 2) ----------------------
# Blanco model effective parameters, For calibration curve calculations (95% CI)
SPEC = dict(
    serial="ATT-11-16-CA85 S/N #02721",
    model="THz Tunable Precision Attenuator ATT-11-16-CA85",
    aperture_mm=85.0,              # clear aperture (holder diameter)
    P_nominal_um=16.0,
    D_phys_um=11.0,                # nominal wire diameter
    P_geom_um=(15.93, 5.03),       # measured via micrograph, quality class "A"
    P_eff_um=(15.50, 0.05),        # 95 % CI
    D_eff_um=(5.67, 0.08),         # 95 % CI, drift-filtered ensemble mean
    loss_db=(0.255, 0.035),        # dB/THz^gamma, 95 % CI  (= alpha 0.0587 Np * 4.343)
    gamma=1.58,
    theta_offset_deg=(-0.45, 0.08),  # 95 % CI, rotator backlash compensation
    band_thz=(0.05, 2.0),          # calibration frequency range
    model_rmse_db=1.39,            # model approximation error, up to -45 dB
    scale_division_deg=1.0,        # цена деления шкалы (подтверждено владельцем)
)

K95 = 1.959963985
# Наши фиты: какие ключи _calib_raw.json соответствуют паспортным параметрам
RAW_KEYS = {"D_eff_um": "D_um", "loss_db": "loss_factor",
            "theta_offset_deg": "angle_offset"}


def session_spread(raw, key):
    """1 sigma по разбросу между сеансами (или None, если данных нет)."""
    sess = [k for k in raw if not k.startswith("_")]
    if len(sess) < 2 or key not in raw[sess[0]]:
        return None
    v = np.array([raw[s][key]["v"] for s in sess], dtype=float)
    return float(v.std(ddof=1))


def conservative(spec_val, spec_ci95, raw, raw_key):
    """Значение из паспорта + max(паспортная сигма, сигма по сеансам)."""
    sigma_spec = spec_ci95 / K95
    sigma_sess = session_spread(raw, raw_key) if raw_key else None
    if sigma_sess is None:
        return Param(spec_val, sigma_spec, "spec_95ci")
    if sigma_sess > sigma_spec:
        return Param(spec_val, sigma_sess, "session_spread")
    return Param(spec_val, sigma_spec, "spec_95ci")


def build(raw_path=None) -> Passport:
    raw_path = Path(raw_path or PASSPORTS / "_calib_raw.json")
    raw = json.loads(raw_path.read_text(encoding="utf-8")) if raw_path.exists() else {}
    sessions = [k for k in raw if not k.startswith("_")]

    p = Passport(
        serial=SPEC["serial"], model=SPEC["model"], aperture_mm=SPEC["aperture_mm"],
        P_um=Param(SPEC["P_eff_um"][0], SPEC["P_eff_um"][1] / K95, "spec_95ci"),
        D_eff_um=conservative(*SPEC["D_eff_um"], raw, RAW_KEYS["D_eff_um"]),
        D_phys_um=SPEC["D_phys_um"],
        loss_factor=conservative(*SPEC["loss_db"], raw, RAW_KEYS["loss_db"]),
        gamma=Param(SPEC["gamma"], 0.0, "spec_fixed"), gamma_range=(1.0, 2.5),
        angle_offset_deg=conservative(*SPEC["theta_offset_deg"], raw,
                                      RAW_KEYS["theta_offset_deg"]),
        tau_ps=Param(0.0, 0.0, "fixed"), tau_par_ps=Param(0.0, 0.0, "fixed"),
        band_thz=SPEC["band_thz"],
        fit_scheme="S0",           # паспортная формула: один вращаемый WGP + проекция
        fit_detector="coherent",   # Scenario A (PCA Menlo Tera K8)
        datasets=tuple(sessions),
        scale=ScaleSpec(division_deg=SPEC["scale_division_deg"]),
        film=film_for_aperture(SPEC["aperture_mm"]),
        model_rmse_db=SPEC["model_rmse_db"],
        model_rmse_note=("паспортная RMSE аппроксимации, динамический диапазон до -45 дБ; "
                         "доминируют неучтённые переотражения между решётками и "
                         "дифракционное рассеяние при theta > 50 град"),
        calibration_note=(
            f"Source: official product passport {SPEC['serial']} "
            f"(passport.tex), Menlo Tera K8 THz-TDS, horizontal linear "
            f"polarization. Device = TWO independent wire-grid polarizers; "
            f"the first is aligned with the beam polarization (theta1=0), the second "
            f"(analyzer) rotates. Attenuation is set by the relative angle theta=theta2-theta1. "
            f"Nominal P/D = {SPEC['P_nominal_um']}/{SPEC['D_phys_um']} um; "
            f"micrograph gives P_geom = {SPEC['P_geom_um'][0]}+-{SPEC['P_geom_um'][1]} um. "
            f"D_eff = {SPEC['D_eff_um'][0]} um is an equivalent FLAT strip "
            f"(a round wire D=11 acts roughly as a strip of width D/2), "
            f"not the geometric diameter. "
            + (f"Sessions used for reproducibility estimate: {', '.join(sessions)}."
               if sessions else "")),
    )
    for w in p.validate():
        print(f"  warning: {w}")
    return p


def main():
    p = build()
    out = PASSPORTS / "ATT-11-16-CA85_02721.json"
    p.save(out)
    print(f"=== passport {p.serial} ===")
    print(f"  aperture   {p.aperture_mm:g} mm -> scale {p.scale.division_deg:g} deg, "
          f"sigma(angle) {p.scale.sigma_total_deg():.2f} deg")
    for name, prm, unit in (("P_eff", p.P_um, "um"), ("D_eff", p.D_eff_um, "um"),
                            ("loss", p.loss_factor, f"dB/THz^{p.gamma.value:g}"),
                            ("offset", p.angle_offset_deg, "deg")):
        print(f"  {name:<10} {prm.value:+8.3f} +- {prm.sigma:.3f} {unit:<18} [{prm.source}]")
    print(f"  D_eff/D_phys = {p.D_eff_um.value / p.D_phys_um:.3f}")
    print(f"  model RMSE   {p.model_rmse_db:.2f} dB")
    print(f"  band         {p.band_thz[0]}-{p.band_thz[1]} THz, scheme {p.fit_scheme}, "
          f"detector {p.fit_detector}")
    print(f"  anomaly nu   {p.nu_anomaly_thz:.1f} THz, green zone up to {p.zone_green_max_thz:.2f} THz")
    print(f"  films        {p.film.name} (optional in scheme S3)")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
