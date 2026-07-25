"""Generate the Monday PureWave experiment-protocol PDF (Russian, Cyrillic-safe)."""
import sys
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                ListFlowable, ListItem, HRFlowable)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

OUT = Path(__file__).resolve().parents[1] / "EXPERIMENT_PROTOCOL_purewave.pdf"

# --- Cyrillic fonts (Windows Arial) ---
FONTS = Path("C:/Windows/Fonts")
pdfmetrics.registerFont(TTFont("Body", str(FONTS / "arial.ttf")))
pdfmetrics.registerFont(TTFont("Body-Bold", str(FONTS / "arialbd.ttf")))

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontName="Body-Bold", fontSize=15, spaceAfter=6, textColor=colors.HexColor("#1a3b5d"))
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontName="Body-Bold", fontSize=12, spaceBefore=8, spaceAfter=4, textColor=colors.HexColor("#245080"))
P  = ParagraphStyle("P", parent=ss["BodyText"], fontName="Body", fontSize=9.5, leading=13, spaceAfter=4)
PS = ParagraphStyle("PS", parent=P, fontSize=8.5, textColor=colors.HexColor("#444444"))
LI = ParagraphStyle("LI", parent=P, spaceAfter=2)

story = []


def h1(t): story.append(Paragraph(t, H1))
def h2(t): story.append(Paragraph(t, H2))
def p(t):  story.append(Paragraph(t, P))
def ps(t): story.append(Paragraph(t, PS))
def sp(h=4): story.append(Spacer(1, h))
def rule(): story.append(HRFlowable(width="100%", thickness=0.6, color=colors.HexColor("#cccccc"), spaceBefore=4, spaceAfter=6))
def bullets(items, style=LI):
    story.append(ListFlowable([ListItem(Paragraph(i, style), leftIndent=10) for i in items],
                              bulletType="bullet", start="•", leftIndent=12))


h1("Протокол эксперимента — угловые THz-TDS измерения PureWave")
ps("Автосгенерировано исследовательским циклом THz-WGP · образец PureWave (P=25.5 µm, D≈10.6 µm, D/P≈0.42) · "
   "цель понедельника: получить первый образец с честной статистикой на точку.")
rule()

h2("1. Зачем это (научная цель)")
p("Текущая модель THz-WGP имеет несколько взаимно-<b>вырожденных</b> частотно-зависимых каналов: "
  "эффективный диаметр D_eff, утечка параллельной компоненты η(ν), не-Друде показатель потерь γ. "
  "Из-за вырождения (a) амплитуда утечки η не экстраполируется с ярких углов на тень; "
  "(b) γ не идентифицируется как физическая константа. Разорвать это можно только "
  "<b>внешними данными</b>: многоповторной статистикой (веса/ошибки) и прямыми доп. измерениями.")
bullets([
  "<b>Многоповторность</b> → σ на каждую точку (θ,ν) → корректный взвешенный χ² и bootstrap-ошибки на η, τ_leak, D_eff.",
  "<b>Тонкая выборка у скрещённого + за скрещённым</b> → локализация точки экстинкции и рычаг на t_par "
  "(у Specac истинный минимум был на ~84°, не 90°).",
  "<b>Доп. эксперименты</b> (если позволяет установка) → прямой тест «γ = диффузное рассеяние на намотке» и "
  "разделение утечки образца от кросс-поляризации детектора.",
])

h2("2. Что понадобится / чек-лист перед стартом")
bullets([
  "Опорный трейс (зеркало/открытый пучок) — сохранить.",
  "Проверить динамический диапазон: трейс у скрещённого (θ≈84–90°) должен быть <b>&gt;10×</b> над шумом; "
  "если нет — усреднять больше повторов в этой зоне.",
  "Записать влажность и температуру в лаборатории на старте и в конце (для контроля дрейфа/линий воды).",
  "Формат файла: два столбца ASCII (время_пс, поле). 3-й столбец (позиция каретки) допустим — игнорируется.",
])

h2("3. Бюджет времени и оптимизация (ВАЖНО)")
p("Один трейс = <b>6 мин 15 с</b> (окно скана 130 пс). Равномерный план из 5 повторов (210 трейсов) = "
  "<b>~22 часа ≈ 3 дня</b> — в один день НЕ влезает. Поэтому план оптимизирован под <b>один рабочий день</b>: "
  "повторы сконцентрированы у скрещённого (там вся физика), фон — периодический.")
bullets([
  "<b>~54 сигнала + ~10 фонов ≈ 64 трейса → ~6.6 ч</b> @6:15 (влезает в день, плотно).",
  "<b>Рычаг №1 — укоротить окно скана до ~65 пс:</b> пик импульса на 25 пс, за ±50 пс лишь 6% энергии "
  "(хвост водяных линий, который мы всё равно маскируем). Разрешение 0.015 ТГц — с запасом для 0.2–1.5 ТГц; "
  "эхо от Δz (~1.3 пс) сохраняется. Ход подвижки ∝ окну → <b>~3:10/трейс → ~3.3 ч</b>; тогда подними повторы "
  "у скрещённого до 6–8.",
  "<b>Рычаг №2 — фон не на каждый трейс:</b> мерить каждые ~2 угла / ~30 мин, копировать в соседние _bg.",
  "<b>Рычаг №3:</b> если 6:15 включает сильное усреднение — снизь его и добери статистику повторами "
  "(та же точность + получаешь σ на точку).",
])
sp(2)
data = [
  ["Зона", "Углы (град)", "Повт."],
  ["Яркая (грубо)", "0, 20, 40, 55, 65, 72", "1"],
  ["У скрещённого (ГЛАВНОЕ, 2°)", "80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100", "4 (6–8 если окно короче)"],
  ["За скрещённым", "105, 110", "2"],
]
t = Table(data, colWidths=[50*mm, 82*mm, 28*mm])
t.setStyle(TableStyle([
  ("FONT", (0,0), (-1,-1), "Body", 8.5),
  ("FONT", (0,0), (-1,0), "Body-Bold", 8.5),
  ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#dbe6f0")),
  ("GRID", (0,0), (-1,-1), 0.4, colors.HexColor("#aaaaaa")),
  ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
  ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f4f7fa")]),
]))
story.append(t); sp(6)
p("Порядок (устойчив к дрейфу): сначала весь набор углов для rep1 (по возрастанию), затем rep2 "
  "(по убыванию) и т.д. — так повторы каждого угла разнесены по времени. Точный список — "
  "<font name='Body-Bold'>data_pool/purewave/MEASUREMENT_ORDER.md</font>.")

h2("4. Протокол на каждый угол")
story.append(ListFlowable([
  ListItem(Paragraph("Повернуть WGP на угол.", LI)),
  ListItem(Paragraph("Записать СИГНАЛ → <font name='Body-Bold'>purewave_&lt;угол&gt;deg_rep&lt;N&gt;_sig.txt</font>.", LI)),
  ListItem(Paragraph("Измерять ФОН (пустой пучок) минимум каждые 2 угла → <font name='Body-Bold'>..._bg.txt</font>. "
                     "Если держишь один фон на блок — скопируй его в _bg каждого повтора этого блока.", LI)),
  ListItem(Paragraph("У скрещённого (80–100°) при слабом сигнале — усредняй больше (добавь rep6, rep7… тем же шаблоном имени).", LI)),
], bulletType="1", leftIndent=14))

h2("5. Дополнительные измерения — если позволяет установка (высокая ценность)")
p("<b>5a. Прямой тест «γ ↔ неоднородность намотки».</b> Поставь приёмник на несколько градусов "
  "ОТ зеркального направления и померь <b>диффузно-рассеянную</b> мощность vs частота (можно при паре "
  "фиксированных углов WGP, напр. 0° и 45°). Если она растёт с частотой как ν^γ с тем же γ, что в фите — "
  "это <b>прямое доказательство</b>, что γ = рассеяние на неоднородности намотки (а не фиктивный параметр). "
  "Именно это разрывает вырождение γ↔утечка.")
p("<b>5b. Разделение утечки ОБРАЗЦА vs кросс-поляризации ДЕТЕКТОРА.</b> Повтори короткий угловой скан "
  "(хотя бы 0/45/84/90°) с эталонным поляризатором перед детектором (или повернув детекторный канал). "
  "Утечка образца зависит от угла вращения WGP; детекторная — нет. Разность конфигураций их разделит "
  "(сейчас в анализе они разделены лишь частично, corr 0.73).")

h2("6. Именование файлов (стандарт)")
p("<font name='Body-Bold'>purewave_&lt;угол&gt;deg_rep&lt;N&gt;_&lt;sig|bg&gt;.txt</font>  — угол целый (0…110), "
  "N — номер повтора. Пустые файлы уже созданы в data_pool/purewave/ — просто записывай в них. "
  "Загрузчик подхватит датасет автоматически как «purewave».")

h2("7. Задание на дополнительные микрофотографии")
bullets([
  "<b>356att (D11/P16) — КРИТИЧНО, сейчас заблокировано:</b> нужна <b>калибровочная граткула на ТОМ ЖЕ "
  "увеличении</b>, что снимки этой сетки (перенос калибровки PureWave не сработал — период вышел 40–94 µm "
  "вместо 16). Без своей граткулы геометрия 356att не извлекается → нет 4-й точки закона D_eff.",
  "<b>Диаметр всех образцов:</b> снимок на <b>максимальном</b> увеличении по <b>резким</b> проволокам "
  "(блик занижает D ~вдвое; истинный D = ширина следа на 10% высоты).",
  "<b>PureWave highmag</b> для подтверждения D (сейчас 10.6 µm — по одному методу).",
  "Имена: &lt;sampleID&gt;__cal__&lt;low|mid|high&gt;mag__div0p01mm.bmp и "
  "&lt;sampleID&gt;__wgp__&lt;mag&gt;mag[_NN].bmp; граткула и WGP на ОДНОМ увеличении.",
])

h2("8. Что это даст (критерий успеха понедельника)")
bullets([
  "Первый образец с σ на точку → взвешенный фит и реальные ошибки на η, τ_leak, D_eff.",
  "Проверка: сужаются ли ошибки у скрещённого (там живёт вся физика утечки).",
  "Если сделаны 5a/5b — прямой ответ на вопрос про γ и разделение образец/детектор.",
  "PureWave (D/P=0.42) станет 4-й точкой закона D_eff (после test_grid 0.46, specac 0.56, 356att 0.71).",
])

ps("Файлы-дорога: data_pool/purewave/*.txt (210 шт.) + MEASUREMENT_ORDER.md. "
   "Аналитический скрипт взвешенного фита будет готов к приёму этих данных.")


def build():
    doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=16*mm, rightMargin=16*mm,
                            topMargin=14*mm, bottomMargin=14*mm, title="PureWave experiment protocol")
    doc.build(story)
    print(f"PDF -> {OUT}  ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    build()
