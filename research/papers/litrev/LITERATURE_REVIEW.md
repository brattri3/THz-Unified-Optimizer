# Литературный обзор: аппроксимация THz-измерений проволочных поляризаторов

**Составлен**: 2026-07-23 | **Метод**: цитатные цепочки через OpenAlex API (references назад + citing works вперёд) от 3 затравочных статей + тематические поиски по ключевым словам.
**Сырые данные**: `*_refs.json`, `*_citedby.json` в этой же папке.

> Semantic Scholar API отдавал 429 (троттлинг анонимного пула) — весь обзор построен на OpenAlex (без ключа, надёжен). Connected Papers открытого API не имеет.

---

## 0. Три затравочные статьи (якоря темы)

| # | Работа | Роль в проекте | OpenAlex |
|---|---|---|---|
| **S1** | **Blanco A., Fonti S., Piacente A. (1986)**, *Transmission coefficients of free-standing wire grids in the far infrared: a theoretical approach for easy computation*, **Infrared Physics 26(6):357–363**, DOI [10.1016/0020-0891(86)90058-8](https://doi.org/10.1016/0020-0891(86)90058-8) | **Аналитическое ядро модели** (эквивалентные схемы, t⊥/t∥); границы: $d/\lambda<0.2$, $d/p<0.1$ (приемлемо до 0.5) | индексируется в OpenAlex и Crossref (запись исправлена L 2026-07-29, см. `answers/D2_P1_blanco1986.md`) |
| **S2** | **Manabe & Murk (2005)**, *Transmission and Reflection Characteristics of Slightly Irregular Wire-Grids with Finite Conductivity…*, IEEE TAP 53(1):250 | **Конечная проводимость + теория возмущений для нерегулярности/шероховатости** (гипотеза диффузного рассеяния) | W2138091822 |
| **S3** | **Castro-Camus (2011/2012)**, *Polarization-Resolved Terahertz Time-Domain Spectroscopy*, JIMTW 33:418 | **Поляриметрия THz-TDS + многоконтактные PCA-детекторы** (гипотеза кросс-поляризации детектора) | W2100964144 |

Затравки покрывают ровно три столпа проекта: **(теория решётки) + (физика потерь) + (аппаратная поляриметрия)**.

---

## 1. Теоретическая ветвь: электродинамика проволочных решёток

Родословная модели (из references S2 и форвардных цитирований):

- **Wait, J. R. (1955)** — *Reflection at arbitrary incidence from a parallel wire grid* (c=73). Классика: импедансная модель решётки при произвольном угле.
- **Larsen (1962)** — *A Survey of the Theory of Wire Grids* (c=123). Обзорная база эквивалентных схем — фундамент формул Бланко.
- **Yasumoto (1999)** — *Efficient calculation of lattice sums for free-space periodic Green's function* (c=140) и **Electromagnetic Scattering from Periodic Arrays of Two Circular Cylinders (2000)** (c=62). Точный lattice-sum метод, который Manabe&Murk расширили на произвольный угол — **точная альтернатива аналитике Бланко** для валидации.
- Современные продолжения (citing S2):
  - **Effects of Random Positioning Errors Upon EM Characteristics of a Wire Grid (2011, c=10)** — прямое развитие темы нерегулярности шага (релевантно вашей гипотезе «эффективного периода»).
  - **Wire grid embedded in SNG/SZ/chiral/bi-isotropic media (2014–2019)** — обобщения импедансной границы (потенциально для подложечных поляризаторов).
  - **Analytical Modeling of Microwave Transmission Through Metasurfaces… Wire Grids (2024, c=1)** — свежая аналитика метаповерхностей.
  - **Wire-Grid and Sparse MoM Antennas: Past/Present/Future (2023, c=32)** — MoM-обзор, мост к численным методам.

**Вывод**: точный численный ориентир для проверки «эффективного диаметра/периода» — это **lattice-sum метод Yasumoto/Manabe**, а не только аналитика Бланко. Кандидат №1 для сверки D_eff.

---

## 2. Физика потерь: конечная проводимость, скин-эффект, шероховатость

- **Manabe & Murk (2005, S2)** — перенос нерегулярности шага в возмущение к регулярной теории; показано, что нерегулярность влияет сильно на **TM (E∥проводам)** и почти не влияет на TE. ⚠️ Это прямо пересекается с вашей гипотезой τ_par (анизотропия при θ→90°).
- **Laman & Grischkowsky (2008)**, *Terahertz conductivity of thin metal films*, APL 93:051105 — отклонение ТГц-проводимости от Друде из-за границ зёрен/шероховатости; обоснование дробного показателя γ (между скин-эффектом ν^0.5 и Рэлеем ν^4). Цитируется в `scientific_foundation.md`.
- Ветвь свободновисящих вольфрамовых решёток (references S2/S3) — **прямые аналоги вашего объекта**:
  - *Free-standing grids wound from 5 μm diameter wire for far-IR spectroscopy (1979, c=37)*
  - *Performance of free-standing grids wound from 10-μm tungsten wire at submm wavelengths: computation and measurement (1981, c=16)*
  - *Fabrication and Characterization of Large Free-Standing Polarizer Grids for mm-Waves (1999, c=16)*
  - *Fabrication of WGP and dependence of submm optical performance on pitch uniformity (1990, c=19)* — **эмпирика влияния неоднородности шага** на характеристики.

**Вывод**: канал «омические Друде + диффузное рассеяние ν^γ» имеет прочную опору (Laman-Grischkowsky + Manabe). Работы 1979–1999 по вольфрамовым решёткам — прямые экспериментальные предшественники, обязательны к цитированию в вашей статье.

---

## 3. Экспериментальные THz-WGP (для сравнения регимов)

Из тематического поиска — в основном **фабрикованные плёночные/решёточные** поляризаторы (иной режим, чем свободновисящая вольфрамовая сетка), но важны для контекста «extinction ratio»:

| Год | Работа | c | Чем полезна |
|---|---|---|---|
| 2009 | Terahertz wire-grid polarizers with micrometer-pitch Al gratings | 206 | Эталон характеризации THz-WGP |
| 2012 | Extremely high extinction ratio THz broadband polarizer, bilayer subwavelength | 69 | Предел динамического диапазона (ваши −40 дБ) |
| 2014 | High extinction ratio thin-film THz polarizer, tunable bilayer | 55 | Двухслойные структуры |
| 2012 | Extraordinary optical transmission/extinction in THz WGP | 46 | Резонансные эффекты |
| 2017 | Fabrication of WGP from UV to THz | 31 | Кросс-диапазонная методология |
| 2011 | WGP sheet by nanoimprint | 72 | Технологии массового производства |
| **2018** | **Fabry-Pérot interferometer with WGP as beamsplitters at THz** | 9 | **Прямо к вашей будущей гипотезе Фабри-Перо между решётками** |
| 2010 | Double Wire-Grid THz Polarizer on Low-Loss Polymer | 18 | Двойная решётка (конфигурация A) |

**Вывод**: ваша ниша — **не изготовление, а метрология/обратная задача для свободновисящей вольфрамовой сетки в плотном режиме D/P>0.5**. В этом кластере таких работ почти нет → подтверждает заявленную научную новизну.

---

## 4. Обратная задача THz-TDS (методология извлечения параметров)

Из тематического поиска — прямая опора для вашего оптимизационного подхода:

- **Duvillaret et al. (1996)** — *A reliable method for extraction of material parameters in THz-TDS* (c=970). **Канонический** алгоритм инверсии; уже в `scientific_foundation.md`.
- **Fixed-point iteration extraction (2005, c=54)** — устойчивость инверсии.
- **Self-calibrating technique for THz-TDS extraction (2011, c=17)** — самокалибровка (перекликается с вашей идеей real-time reference tracking).
- **ANN for material parameter extraction in THz-TDS (2022, c=31)** — ML-инверсия, возможное расширение вашего оптимизатора.
- **Transmission vs Reflection extraction comparative (2018, c=24)**.

**Вывод**: позиционируйте свой 2D-фитинг Бланко как расширение линии Duvillaret на **геометрические** (не только материальные) параметры решётки. ANN-2022 — потенциальный next step против «flat valley» вырождения.

---

## 5. Аппаратная поляриметрия и кросс-поляризация детектора

Прямо под вашу гипотезу об аппаратной функции (ε_cross, θ_offset):

- **Multi-contact photoconductive receivers** (references S3): *Polarization-sensitive THz detection by multicontact PCA (2005, c=144)*, *3-contact PCA (2007, c=64)*, *4-contact PCA (2006, c=75)*.
- **A matter of symmetry: THz polarization detection of multi-contact PCA via response matrix analysis (2015, c=16)** — формализм матрицы отклика детектора; **напрямую моделирует паразитную кросс-поляризацию**.
- **A polarization-sensitive 4-contact detector for THz-TDS (2014, c=23)**.
- **Constraints on Jones transmission matrices from time-reversal & spatial symmetries (2014, c=73)** — строгие ограничения на матрицы Джонса (проверка корректности вашей матричной модели).
- **THz-TDP spinning E-O sampling: precision & calibration (2020, c=34)** и **highly precise/accurate THz polarization via EO with modulation (2014, c=51)** — метрология точности поляризации, эталон для вашего заявления о «пределе точности прибора».

**Вывод**: гипотеза «кросс-поляризация детектора даёт ложный сигнал на 90°» имеет прямую опору в response-matrix формализме (2015) и multi-contact PCA. **Обязательно к цитированию** — это ядро вашей аппаратной поправки.

---

## 6. Карта «гипотеза проекта → литературная опора»

| Гипотеза проекта | Ключевая литература | Статус опоры |
|---|---|---|
| Аналитика Бланко (t⊥/t∥) | Blanco 1986, Larsen 1962, Wait 1955 | прочная |
| Эффективный D_eff / P_eff (D/P>0.5) | Yasumoto/Manabe lattice-sum; positioning errors 2011 | **пробел — ниша новизны** |
| Импеданс Друде (вольфрам) | Laman & Grischkowsky 2008 | прочная |
| Диффузное рассеяние ν^γ | Manabe&Murk 2005; Laman 2008 | прочная |
| Фазовая анизотропия τ_par (TE/TM) | Manabe&Murk (TM≫TE эффект нерегулярности) | косвенная — требует усиления |
| Кросс-поляризация детектора ε_cross | multi-contact PCA; response-matrix 2015; Castro-Camus 2011 | прочная |
| Компенсация дрейфа лазера | Naftaly & Dudley 2009; self-calibrating 2011 | средняя |
| Фабри-Перо между решётками (future) | FP-interferometer WGP beamsplitters 2018 | готовая опора |
| Дифракция на апертуре (future) | — | пробел |

---

## 7. Рекомендованные к прочтению (приоритет)

1. **Manabe & Murk 2005** (уже есть PDF) — перечитать раздел про TM/TE асимметрию нерегулярности → усилить τ_par.
2. **Response-matrix analysis of multi-contact PCA (2015)** — формальная база ε_cross.
3. **Yasumoto lattice-sum (1999/2000)** — точный численный ориентир для валидации D_eff.
4. **Laman & Grischkowsky 2008** — обоснование дробного γ.
5. **FP-interferometer WGP beamsplitters (2018)** — задел под гипотезу Фабри-Перо.
6. **ANN extraction 2022** — против вырождения параметров.

---

## 8. Пробелы и возможности (научная новизна)

- **Плотный режим D/P>0.5 свободновисящих вольфрамовых решёток** методом THz-TDS обратной задачи — практически не покрыт литературой. Это ваша главная ниша.
- **Совместное разделение** физики образца (Бланко+Друде+рассеяние) и аппаратной физики (ε_cross детектора) в единой 2D-оптимизации — методологически ново.
- **Эмпирический закон D_eff/D_phys = 1 − 0.85·(D/P)** — нет прямых аналогов; требует сверки с точным lattice-sum расчётом для публикабельности.
