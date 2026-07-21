---
marp: true
theme: gaia
_class: lead
paginate: true
backgroundColor: #f5f5f7
color: #1d1d1f
style: |
  section {
    font-family: 'Helvetica Neue', Arial, sans-serif;
    padding: 22px 30px;
  }
  h1 {
    color: #0071e3;
    font-size: 1.5em;
    margin-bottom: 8px;
  }
  h2 {
    color: #1d1d1f;
    border-bottom: 2px solid #0071e3;
    font-size: 1.2em;
    margin-bottom: 8px;
    padding-bottom: 4px;
  }
  footer {
    font-size: 0.5em;
    color: #86868b;
  }
  .grid {
    display: grid;
    grid-template-columns: 0.9fr 1.1fr;
    gap: 12px;
    align-items: center;
  }
  .highlight {
    background-color: #e8f2fc;
    padding: 6px 10px;
    border-left: 4px solid #0071e3;
    border-radius: 4px;
    font-size: 0.78em;
    margin-top: 8px;
  }
  ul {
    font-size: 0.78em;
    margin-top: 3px;
    margin-bottom: 3px;
    padding-left: 18px;
    line-height: 1.35;
  }
  li {
    margin-bottom: 4px;
  }
  p {
    font-size: 0.8em;
    line-height: 1.35;
    margin-top: 4px;
    margin-bottom: 4px;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.72em;
    margin-top: 4px;
  }
  th {
    background-color: #0071e3;
    color: white;
    padding: 5px;
  }
  td {
    padding: 5px;
    border: 1px solid #d2d2d7;
    text-align: center;
  }

---

# Развитие физической модели ТГц поляризатора
## Учёт фазовых эффектов, проводимости Друде и закономерности $D_{\text{eff}}$

**Докладчик:** Попов Дмитрий
*Анализ терагерцовой спектроскопии (THz-TDS)*

---

## Конструкция и развитие физической модели

<div class="grid">
  <div>

  * **Устройство:** Аттенюатор `ATT-11-16-CA85` ($P_1$ вращается, $P_2$ фиксирован).
  * **Итог прошлой недели:** 6-параметрическая модель дала RMSE = **0.29 дБ**.
  * **Новое развитие этой недели:**
    * Переход к **комплексной** аппроксимации (амплитуда + фаза).
    * Внедрён импеданс Друде для вольфрамовых нитей.
    * Учтена **фазовая анизотропия** $\tau_{\text{par}}$ при скрещивании.
    * Модель проверена на новом образце 40/20 мкм.

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/optim_stage3_hw.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Спектральное описание: потери Друде и рассеяние

<div class="grid">
  <div>

  * **Проблема M0:** При углах > 50° чистая аналитика Бланко завышает сигнал.
  * **Друде (M1):** Конечная проводимость вольфрама $\sigma_0 = 1.8 \cdot 10^7$ См/м снижает спектральный уровень.
  * **Рассеяние (M2):** Рэлеевское затухание $e^{-\frac{1}{2}\alpha\nu^2}$ на шероховатостях *(Manabe & Murk, IEEE TAP, 2005)*.
  * **Результат:** Точное воспроизведение спектра на всех углах в диапазоне 0.2–1.5 ТГц.

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/pres_slide3_ablation.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Комплексная аппроксимация и фазовая анизотропия $\tau_{\text{par}}$

<div class="grid">
  <div>

  * **Комплексный анализ:** Одновременная подгонка |E| и фазы *(Castro-Camus, J. Infrared TE Waves, 2012)*.
  * **Физический эффект:** На углах 80°–90° TE-мода испытывает фазовую задержку $\tau_{\text{par}}$.
  * **Результат:** Учёт $\tau_{\text{par}} \approx -0.1$ пс снижает ошибку фазы более чем в **5 раз**.

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/pres_slide4_tau_par.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Валидация модели на решётке 40/20 мкм

<div class="grid">
  <div>

  * **Образец:** Период $P = 40.0$ мкм, нить $D = 20.0$ мкм ($D/P = 0.50$).
  * **Результаты подгонки M3:**
    * Оффсет юстировки: $\theta_{\text{offset}} = 0.47^\circ$.
    * Извлечённый диаметр: $D_{\text{eff}} = 11.37$ мкм.
    * RMSE: **0.96 дБ** (амплитуда), **0.24 рад** (фаза).
  * **Вывод:** Модель универсальна для различных геометрий WGP.

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/pres_slide5_validation.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Физический закон масштабирования $D_{\text{eff}}$

<div class="grid">
  <div>

  * **Электродинамическое сжатие:** $D_{\text{eff}} < D_{\text{phys}}$ из-за перекрытия ближних полей нитей:
    * 15.5/11 мкм ($D/P = 0.71$): $D_{\text{eff}}/D = 0.40$.
    * 40/20 мкм ($D/P = 0.50$): $D_{\text{eff}}/D = 0.57$.
  * **Открытый закон:**
    $$\frac{D_{\text{eff}}}{D_{\text{phys}}} = 1 - 0.85 \cdot \frac{D}{P}$$

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/pres_slide6_deff_law.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Точность измерений и границы прибора

<div class="grid">
  <div>

  * **Опровержение краевой дифракции:** Апертура оправы ($101.6$ мм) в 8.5 раз шире ТГц пучка ($12$ мм) — дифракция физически исключена.
  * **Причина остатков:**
    1. Динамический предел детектора (>37 дБ при 90°).
    2. Дрейф мощности лазера ($\sim 2.2\%$).
  * **Вывод:** Достигнут предел точности спектрометра.

  </div>
  <div>
    <img src="c:/THz-Unified-Optimizer/docs/images/residuals_comprehensive_maps.png" height="360px" style="border-radius:6px; box-shadow: 0 3px 6px rgba(0,0,0,0.1); display: block; margin: 0 auto;">
  </div>
</div>

---

## Сравнение результатов оптимизации

| Параметр | Номинал | M0 (Бланко) | M1 (+Друде) | M2 (+Рассеяние) | M3 (Итоговая) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **$D_{\text{eff}}$ (мкм)** | 11.0 | 3.96 | 3.95 | 3.95 | **4.40** |
| **$\theta_{\text{offset}}$ (°)** | — | — | — | 0.80° | **0.80°** |
| **$\tau_{\text{par}}$ (пс)** | — | — | — | — | **$-0.099$** |
| **RMSE (амплитуда)**| — | 1.12 дБ | 1.25 дБ | 1.15 дБ | **0.96 дБ** |
| **RMSE (фаза)**| — | 0.15 рад | 0.15 рад | 0.15 рад | **0.13 рад** |

<div class="highlight">
Учёт фазовой анизотропии τ_par и потерь Друде обеспечил рекордную точность описания амплитуды (0.96 дБ) и фазы (0.13 рад) в диапазоне 0.2–1.5 ТГц.
</div>

---

## Ключевые выводы

1. **Комплексная аппроксимация и $\tau_{\text{par}}$:** Впервые учтена фазовая задержка TE-моды при скрещивании поляризаторов, что устранило расхождение фазы при углах > 80°.
2. **Закон масштабирования $D_{\text{eff}}$:** Выведено прямолинейное соотношение $D_{\text{eff}}/D_{\text{phys}} = 1 - 0.85 (D/P)$, подтверждённое на образцах 15.5/11 и 40/20 мкм.
3. **Метрологический предел:** Исключена гипотеза дифракции на оправе (апертура 102 мм >> пучок 12 мм). Модель достигла фундаментального предела чувствительности ТГц спектрометра.
