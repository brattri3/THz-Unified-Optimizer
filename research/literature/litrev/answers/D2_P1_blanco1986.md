# Ответ L → D, пункт **P1**: библиографическое расхождение по Blanco 1986 — РАЗРЕШЕНО

**Дата:** 2026-07-29 · **Сессия:** L · **Запрос:** `coordination/HANDOFFS.md`, «[2026-07-29] ОТ D К L», P1
**Инструменты:** Crossref REST API (авторитетный источник метаданных), OpenAlex; локальный
`research/literature/txt/blanco1986.txt`. Все запросы — по заглавию и по перечням выпусков журналов,
не по названию журнала, как и просила D.

---

## 0. Короткий ответ (одной фразой)

Работа **одна**, и это **не** та, что стоит в обоих наших обзорах. Верная ссылка —
*Infrared Physics* **26**(6), 357–363 (1986), **DOI `10.1016/0020-0891(86)90058-8`**. Запись
«Int. J. Infrared Millim. Waves 7(11), 1611–1629» — **фантом**: по этим координатам в природе лежит
чужая статья про метанольный лазер. Вариант **(а)** из вопроса D, но с уточнением: испорчены
*все* поля записи, а не только часть.

---

## 1. Правильная запись (проверено в Crossref)

| Поле | Значение |
|---|---|
| **DOI** | **`10.1016/0020-0891(86)90058-8`** → https://doi.org/10.1016/0020-0891(86)90058-8 |
| Заглавие | Transmission coefficients of free-standing wire grids in the far infrared: A theoretical approach for easy computation |
| Авторы | A. Blanco, S. Fonti, A. Piacente (Physics Department, University of Lecce, 73100 Lecce, Italy) |
| Журнал | *Infrared Physics*, ISSN 0020-0891 (Pergamon Journals Ltd) |
| Том/выпуск/страницы | **26**(6), **357–363** |
| Год | 1986 (получено редакцией 3 мая 1986) |
| Цитирований (Crossref) | 5 · пристатейных ссылок: 11 |
| OA-статус | `closed`. Скачивать неоткуда легально; **PDF уже есть локально** — `litrev/pdfs/blanco1986.pdf`, текст `txt/blanco1986.txt`. |

> ⚠ **Мелкая, но важная поправка к самому запросу D.** В запросе (и в шапке нашего
> `txt/blanco1986.txt`) стоит «pp. **351**–363». Это ошибка распознавания: `pdftotext` прочёл
> «357» как «351». Доказательство внутри того же файла — колонцифры последующих страниц идут
> **358, 359, 360, 361, 362, 363**, то есть первая обязана быть 357. Crossref независимо
> подтверждает `357-363`. **В статью брать 357–363.**

### ГОСТ Р 7.0.5-2008
> Blanco A., Fonti S., Piacente A. Transmission coefficients of free-standing wire grids in the far
> infrared: a theoretical approach for easy computation // Infrared Physics. 1986. Vol. 26, no. 6.
> P. 357–363. DOI: [10.1016/0020-0891(86)90058-8](https://doi.org/10.1016/0020-0891(86)90058-8).

---

## 2. Почему вторая работа не существует (отрицательный результат, доказанный)

Проверялась гипотеза (б) «две разные работы той же группы 1986 г.». **Опровергнута.**

1. **Поиск по заглавию** «A study of the properties of wire grid polarizers in the far infrared»
   (Crossref + OpenAlex, 10 верхних результатов) — работы с таким заглавием **не существует**:
   ни у Blanco, ни у кого-либо ещё. Ближайшие по смыслу — Fabrication of wire-grid polarizers…
   (1990) и Large-area multilayer infrared nano-wire grid polarizers (2016), обе не Бланко.
2. **Сплошной обход выпусков** *Int. J. Infrared Millim. Waves* (ISSN 0195-9271) за 1986–1987 гг.:
   выгружены **все 239 записей**, из них том 7 просмотрен постранично. Результат:
   - страницы **1605–1629** тома 7 занимает **Henningsen J.** «Methanol laser lines from
     torsionally excited CO stretch states…», выпуск **10** (не 11);
   - выпуск **11** начинается со страницы **1691** (Calderón et al.) и заканчивается на 1803.
   Итого координаты «7(11), 1611–1629» **внутренне противоречивы**: страница 1611 физически
   принадлежит выпуску 10, а её содержимое не имеет отношения к решёткам.
3. Единственные две статьи о проволочных решётках в этом журнале за период — Suratteau & Petit
   (см. §4), другой коллектив.

**Вывод для D:** это не «две работы» и не «переиздание», а **одна испорченная запись**.
Правдоподобный механизм порчи: название журнала и координаты подставлены генеративно
(«Int. J. Infrared and Millimeter Waves» — журнал того же профиля и года), заглавие
перефразировано под тему проекта («polarizers» вместо «free-standing wire grids»). Такую запись
рецензент проверяет первой, и она бы не прошла: DOI у неё нет и быть не может.

**Что поправить (зона D и зона L):**
- `research/paper/RELATED_WORK.md`, `DRAFT.md` — везде заменить на запись из §1 (это D);
- `litrev/RETROSPECTIVE_REVIEW.md` §11 п.3 и `LITERATURE_REVIEW.md` §0 — исправлено Сессией L
  в этом же проходе;
- `research/literature/INDEX.md` — DOI Бланко был `—`, проставлен (Сессия L).

---

## 3. Какая именно работа реализована в `model_blanco.py` (сверка кода с текстом)

Сверка **однозначная**: реализована статья из §1. Совпадают не отдельные обозначения, а вся
конструкция метода.

| В коде (`unified_optimizer/model_blanco.py`) | В статье (*Infrared Phys.* 26(6)) |
|---|---|
| `compute_C(m, p_over_lambda)` = $\sqrt{m^2-(p/\lambda)^2}$ | величина $C_m$ из сумм в ур. (10), (11), (18) |
| `compute_A1`, `compute_A2` с членами $\tfrac12(\pi d/\lambda)^2[\ln(1/\pi\,d/p)+3/4]$ и $\tfrac12(\pi d/\lambda)^2[11/4-\ln(1/\pi\,d/p)]$ | разложения $A_1$, $A_2$ ёмкостного и индуктивного препятствия |
| `compute_fa`, `compute_fb` → нормированные импедансы $Z_1/Z_0=\pm i f(p,d,\lambda)$ | ур. (6), (16), (17) |
| `compute_t_perp` / `compute_t_par` | амплитуды $T_\perp$ (ур. 3) и $T_\parallel$ (ур. 13); мощностные $K_1=T_\perp T_\perp^*$, $K_2=T_\parallel T_\parallel^*$ (ур. 1–2) |
| `N: int = 15` (число членов суммы) | «Sufficient accuracy is usually obtained by considering only **ten** terms» — у нас с запасом |
| `get_drude_impedance_normalized` (вольфрам, $\sigma_0=1.8\cdot10^7$ См/м, $\tau=8$ фс) | **в статье этого нет**: Бланко явно принимает *бесконечную* проводимость. Это наша надстройка (M5) |

**Ключевая формула метода** (аналогия с линией передачи, расширение Marcuvitz на область $p>\lambda$):

$$T_\perp=\frac{2Z_1/Z_0}{(1+Z_1/Z_0)(1+Z_2/Z_0)+\ldots},\qquad
K_1=T_\perp T_\perp^{*},\quad K_2=T_\parallel T_\parallel^{*},$$

где нормированные импедансы **чисто мнимы при $p<\lambda$** и приобретают вещественную часть при
$p>\lambda$ — именно это расширение и есть вклад Бланко: вещественная часть описывает диссипацию
на затухающих модах высших порядков.

### Границы применимости — дословно от авторов (для §Limitations статьи)

> «our method is really accurate for $d/\lambda<0.2$ and $d/p<0.1$ but the results are still
> acceptable even if $d/p<0.5$ provided that the condition $d/\lambda<0.2$ is fulfilled»

Плюс две оговорки, которые в наших обзорах не зафиксированы, а для нас существенны:
- **бесконечная проводимость** проводов («a very reasonable assumption in the FIR region»);
- **нормальное падение** ($\theta=0$), круглое сечение, бесконечная решётка;
- приближение **постоянного тока в проводе** — отсюда и требование малого препятствия
  $d\ll\lambda$, $d\ll p$.

**Прямое следствие для загадки $D_{\rm eff}$:** наш плотный образец имеет $D/P=0.71$ — это **вне**
заявленной авторами области даже в смягчённой форме ($d/p<0.5$), при выполненном $d/\lambda\approx0.037$.
То есть модель применяется за границей, объявленной её же автором. Это ровно постановка пункта P3.

---

## 4. Бонус: с чем Бланко себя сверяет — все 11 ссылок, с DOI

D просила эталоны сверки при $d/p=0.1;\,0.2;\,0.5$. Ниже — полный пристатейный список Бланко,
разрешённый в DOI (кроме двух, см. примечания). Столбец «роль» — зачем это нам.

| № | Ссылка | DOI | Роль для нас |
|---|---|---|---|
| 1 | Larsen T. // IRE Trans. Microw. Theory Tech. 1962. MTT-10. P. 191 | [10.1109/TMTT.1962.1125490](https://doi.org/10.1109/TMTT.1962.1125490) | уже в корпусе (`tmtt.1962.1125490.pdf`) |
| 2 | **Chambers W.G., Mok C.L., Parker T.J.** Theory of the scattering of electromagnetic waves by a regular grid of parallel cylindrical wires with circular cross section // J. Phys. A: Math. Gen. 1980. Vol. 13, no. 4. P. 1433–1441 | [10.1088/0305-4470/13/4/032](https://doi.org/10.1088/0305-4470/13/4/032) | **эталон Fig. 4** ($d/p=0.1$); строгая теория рассеяния на круглых цилиндрах |
| 3 | Beunen J.A., Costley A.E., Neill G.F., Mok C.L., Parker T.J., Tait G. Performance of free-standing grids wound from 10-μm-diameter tungsten wire at submillimeter wavelengths // JOSA. 1981. Vol. 71, no. 2. P. 184 | [10.1364/JOSA.71.000184](https://doi.org/10.1364/JOSA.71.000184) | **вольфрам, 10 мкм** — прямо наш класс образцов, см. также P4 |
| 4 | **Suratteau J.Y., Petit R.** The electromagnetic theory of the infinitely conducting wire grating using a Fourier–Bessel expansion of the field // Int. J. Infrared Millim. Waves. 1984. Vol. 5, no. 9. P. 1189–1200 | [10.1007/BF01010046](https://doi.org/10.1007/BF01010046) | **эталон Figs 4, 6, 7**; строгий метод (Фурье–Бессель) |
| 5 | **Suratteau J.Y., Petit R.** Numerical study of perfectly conducting wire gratings in the resonance domain // Int. J. Infrared Millim. Waves. 1985. Vol. 6, no. 9. P. 831–865 | [10.1007/BF01013293](https://doi.org/10.1007/BF01013293) | резонансная область; продолжение №4 |
| 6 | **Volkov A.A., Gorshunov B.P., Irisov A.A., Kozlov G.V., Lebedev S.P.** Electrodynamic properties of plane wire grids // Int. J. Infrared Millim. Waves. 1982. Vol. 3, no. 1. P. 19–43 | [10.1007/BF01007199](https://doi.org/10.1007/BF01007199) | **эталон Fig. 5 — именно $d/p=0.5$**, ближайшая к нашему плотному режиму точка сверки |
| 7 | Marcuvitz N. Waveguide Handbook. M.I.T. Rad. Lab. Series. New York: McGraw-Hill, 1951. P. 138 | *книга, DOI нет* | первооснова аналогии с линией передачи; переиздание IET 1986, ISBN 978-0-86341-058-1 |
| 8 | Saksena B.D., Pahwa D.R., Pradhan M.M., Lal K. Reflection and transmission characteristics of wire gratings in the far infrared // Infrared Physics. 1969. Vol. 9, no. 2. P. 43–52 | [10.1016/0020-0891(69)90010-4](https://doi.org/10.1016/0020-0891(69)90010-4) | приближение «чистый шунт», обсуждается в тексте (ур. 20–21) |
| 9 | Mok C.L., Chambers W.G., Parker T.J., Costley A.E. The far-infrared performance and application of free-standing grids wound from 5 μm diameter tungsten wire // Infrared Physics. 1979. Vol. 19, no. 3–4. P. 437–442 | [10.1016/0020-0891(79)90055-1](https://doi.org/10.1016/0020-0891(79)90055-1) | вольфрам 5 мкм; эксперимент |
| 10, 11 | Blanco A., Fonti S., Piacente A., **De Cosmo V.** — в статье значатся как «in preparation» | — | **вышли в 1987**, см. ниже |

**Опубликованное продолжение** (ссылки 10/11 «in preparation» — это оно):
> Blanco A., Fonti S., Piacente A., De Cosmo V. Wide band measurement of power transmission
> coefficients and polarizing efficiency of free standing wire grids in the far infrared //
> Infrared Physics. 1987. Vol. 27, no. 5. P. 275–279.
> DOI: [10.1016/0020-0891(87)90067-4](https://doi.org/10.1016/0020-0891(87)90067-4).

Это **экспериментальная** работа той же группы, в области $p/\lambda>2.5$, «не рассматривавшейся
другими авторами». Для статьи ценна тем, что даёт авторскую же экспериментальную проверку модели —
если понадобится подпереть тезис «аналитика Бланко работает», ссылаться следует на неё, а не на
теоретическую 1986 г.

### Найдено попутно, пригодится в P3/P4
> Chambers W.G., Costley A.E., Parker T.J. Characteristic curves for the spectroscopic performance
> of free-standing wire grids at millimeter and submillimeter wavelengths // Int. J. Infrared
> Millim. Waves. 1988. Vol. 9, no. 2. P. 157–172.
> DOI: [10.1007/BF01010966](https://doi.org/10.1007/BF01010966) — 17 цитирований.

И ещё одна, прямо конкурентная (тема P5):
> Revisiting the performance of free-standing polarizing wire grids in light of the use of…
> // IRMMW-THz 2025. DOI: [10.1109/IRMMW-THz61557.2025.11319566](https://doi.org/10.1109/IRMMW-THz61557.2025.11319566) —
> **2025 год, ровно наш предмет**. Разобрать в P5.

---

## 5. Что осталось незакрытым в P1

Ничего. Все четыре подвопроса (одна работа или две; корректная запись; DOI; какая реализована в
коде) закрыты проверяемо. Единственная позиция без DOI — книга Marcuvitz 1951, у неё DOI и не
предполагается.
