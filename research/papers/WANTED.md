# WANTED — работы, найденные поиском, но недоступные легально из сессии

**Ведёт:** Сессия L · **Обновлено:** 2026-07-29

Здесь то, что найдено и проверено в Crossref/OpenAlex, но чего **нет** в `litrev/pdfs/`, потому что
открытого доступа нет. `metasearch.py` качает только легальный OA (Unpaywall/arXiv/репозитории);
за платную стену он не ходит и ходить не должен.

**Владельцу.** Если у вас есть доступ через организацию — эти файлы можно скачать и положить в
`research/papers/litrev/pdfs/` **под именем из колонки «Имя файла»** (схема `ГОД_Автор_слаг`,
см. `FILENAME_MAP.md`). Текст я извлеку сам следующей командой:

```powershell
# poppler pdftotext, оффлайн; -layout сохраняет колонки
pdftotext -layout "research\papers\litrev\pdfs\<имя>.pdf" "research\papers\txt\<имя>.txt"
```

Приоритет — сверху вниз. Первые три закрывают конкретные открытые пункты запросов D и A.

---

## ★★★ Приоритет 1 — блокируют содержательные ответы

| № | Работа | DOI | Имя файла | Что закроет |
|---|---|---|---|---|
| 1 | **Laman N., Grischkowsky D.** Terahertz conductivity of thin metal films // Appl. Phys. Lett. 2008. Vol. 93, no. 5. Art. 051105 | [10.1063/1.2968308](https://doi.org/10.1063/1.2968308) | `2008_Laman-Grischkowsky_Terahertz-conductivity-thin-metal-films.pdf` | **P2**: формула отклонения $\sigma(\omega)$ от Друде и механизм (границы зёрен). Сейчас тезис статьи о члене $\nu^\gamma$ держится на пересказе. Unpaywall обещает green-копию в SHAREOK, но ссылка **мертва (404)** |
| 2 | **El-Agez T.M., Taya S.A.** An extensive theoretical analysis of the 1:2 ratio rotating polarizer–analyzer Fourier ellipsometer // Physica Scripta. 2011. Vol. 83, no. 2. Art. 025701 | [10.1088/0031-8949/83/02/025701](https://doi.org/10.1088/0031-8949/83/02/025701) | `2011_El-Agez-Taya_1-2-ratio-rotating-polarizer-analyzer-Fourier-ellipsometer.pdf` | **запрос A, п.2**: единственная найденная работа, где высшие гармоники (включая 4ω) разбираются явно. Нужно понять, совпадает ли их механизм с нашим. От этого зависит, можно ли заявлять новизну по 4ω |
| 3 | **Hauge P.S.** Recent developments in instrumentation in ellipsometry // Surface Science. 1980. Vol. 96, no. 1–3. P. 108–140 | [10.1016/0039-6028(80)90297-6](https://doi.org/10.1016/0039-6028(80)90297-6) | `1980_Hauge_Recent-developments-instrumentation-ellipsometry.pdf` | **запрос A, п.1**: канон приёма «лок-ин на 2ω». Нужен, чтобы сослаться на конкретную формулу, а не на факт существования метода |

## ★★ Приоритет 2 — эталоны сверки Бланко (пункт P3, загадка $D_{\rm eff}$)

Это работы, с которыми **сам Бланко** сверяет свою аналитику (его Figs 4–7). Дают независимые
точные расчёты при $d/p = 0.1;\ 0.2;\ 0.5$ — то есть подход к нашему плотному режиму $D/P=0.71$.

| № | Работа | DOI | Имя файла |
|---|---|---|---|
| 4 | **Volkov A.A. et al.** Electrodynamic properties of plane wire grids // Int. J. Infrared Millim. Waves. 1982. Vol. 3, no. 1. P. 19–43 | [10.1007/BF01007199](https://doi.org/10.1007/BF01007199) | `1982_Volkov-etal_Electrodynamic-properties-plane-wire-grids.pdf` |
| 5 | **Chambers W.G., Mok C.L., Parker T.J.** Theory of the scattering of electromagnetic waves by a regular grid of parallel cylindrical wires… // J. Phys. A. 1980. Vol. 13, no. 4. P. 1433–1441 | [10.1088/0305-4470/13/4/032](https://doi.org/10.1088/0305-4470/13/4/032) | `1980_Chambers-etal_Theory-scattering-regular-grid-parallel-cylindrical-wires.pdf` |
| 6 | **Suratteau J.Y., Petit R.** Numerical study of perfectly conducting wire gratings in the resonance domain // Int. J. Infrared Millim. Waves. 1985. Vol. 6, no. 9. P. 831–865 | [10.1007/BF01013293](https://doi.org/10.1007/BF01013293) | `1985_Suratteau-Petit_Numerical-study-perfectly-conducting-wire-gratings.pdf` |
| 7 | **Suratteau J.Y., Petit R.** The electromagnetic theory of the infinitely conducting wire grating using a Fourier–Bessel expansion… // Int. J. Infrared Millim. Waves. 1984. Vol. 5, no. 9. P. 1189–1200 | [10.1007/BF01010046](https://doi.org/10.1007/BF01010046) | `1984_Suratteau-Petit_Electromagnetic-theory-infinitely-conducting-wire-grating.pdf` |

**№ 4 — самый ценный из четырёх:** именно на нём Бланко строит сравнение при $d/p=0.5$, ближайшее
к нашему образцу, что есть в литературе.

## ★★ Приоритет 2 — свободновисящие вольфрамовые решётки (пункт P4, раздел «Образцы»)

| № | Работа | DOI | Имя файла |
|---|---|---|---|
| 8 | **Beunen J.A. et al.** Performance of free-standing grids wound from **10-μm-diameter tungsten wire** at submillimeter wavelengths // JOSA. 1981. Vol. 71, no. 2. P. 184 | [10.1364/JOSA.71.000184](https://doi.org/10.1364/JOSA.71.000184) | `1981_Beunen-etal_Performance-free-standing-grids-10um-tungsten-wire.pdf` |
| 9 | **Blanco A. et al.** Wide band measurement of power transmission coefficients and polarizing efficiency of free standing wire grids // Infrared Physics. 1987. Vol. 27, no. 5. P. 275–279 | [10.1016/0020-0891(87)90067-4](https://doi.org/10.1016/0020-0891(87)90067-4) | `1987_Blanco-etal_Wide-band-measurement-power-transmission-coefficients.pdf` |
| 10 | **Chambers W.G., Costley A.E., Parker T.J.** Characteristic curves for the spectroscopic performance of free-standing wire grids… // Int. J. Infrared Millim. Waves. 1988. Vol. 9, no. 2. P. 157–172 | [10.1007/BF01010966](https://doi.org/10.1007/BF01010966) | `1988_Chambers-etal_Characteristic-curves-spectroscopic-performance-free-standing-wire-grids.pdf` |
| 11 | **Mok C.L. et al.** The far-infrared performance and application of free-standing grids wound from **5 μm tungsten wire** // Infrared Physics. 1979. Vol. 19, no. 3–4. P. 437–442 | [10.1016/0020-0891(79)90055-1](https://doi.org/10.1016/0020-0891(79)90055-1) | `1979_Mok-etal_Far-infrared-performance-free-standing-grids-5um-tungsten.pdf` |

№ 9 — авторская экспериментальная проверка модели Бланко; на неё правильнее ссылаться, чем на
теоретическую 1986 г., когда речь о согласии с экспериментом.

## ★ Приоритет 3 — каноны не-друдевских потерь (P2) и конкуренты (P5)

| № | Работа | DOI | Имя файла |
|---|---|---|---|
| 12 | **Sondheimer E.H.** The mean free path of electrons in metals // Adv. Phys. 1952. Vol. 1. P. 1–42 (переизд. 2001: [10.1080/00018730110102187](https://doi.org/10.1080/00018730110102187)) | [10.1080/00018735200101151](https://doi.org/10.1080/00018735200101151) | `1952_Sondheimer_mean-free-path-electrons-metals.pdf` |
| 13 | **Mayadas A.F., Shatzkes M.** Electrical-resistivity model for polycrystalline films… // Phys. Rev. B. 1970. Vol. 1, no. 4. P. 1382–1389 | [10.1103/PhysRevB.1.1382](https://doi.org/10.1103/PhysRevB.1.1382) | `1970_Mayadas-Shatzkes_Electrical-resistivity-model-polycrystalline-films.pdf` |
| 14 | **Revisiting the performance of free-standing polarizing wire grids…** // IRMMW-THz 2025 | [10.1109/IRMMW-THz61557.2025.11319566](https://doi.org/10.1109/IRMMW-THz61557.2025.11319566) | `2025_IRMMW_Revisiting-performance-free-standing-polarizing-wire-grids.pdf` |
| 15 | **A very simple and accurate way to measure the transmission axis of a linear polarizer** // Meas. Sci. Technol. 2019 | [10.1088/1361-6501/ab3688](https://doi.org/10.1088/1361-6501/ab3688) | `2019_Meas-Sci-Technol_measure-transmission-axis-linear-polarizer.pdf` |

№ 14 — 2025 год, ровно наш предмет (свободновисящие поляризующие решётки), потенциальный конкурент
по приоритету. № 15 — про измерение самого $\theta_0$, смежно с расхождением 9.3σ у Сессии A.

---

## Уже доступно свободно — качать не нужно

- Alhaj Hasan et al. Wire-Grid and Sparse MoM Antennas // Symmetry. 2023. Vol. 15, no. 2. Art. 378 —
  [10.3390/sym15020378](https://doi.org/10.3390/sym15020378), **OA**. Скачаю сам следующим проходом.

## Уже есть в корпусе, но без извлечённого текста

- `1977_Costley-etal_Free-standing-fine-wire-grids-manufacture.pdf` — PDF лежит, `.txt` нет.
  Это **канон изготовления** (пункт P4). Извлечь текст: команда в шапке этого файла.
