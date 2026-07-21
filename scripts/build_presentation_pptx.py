import os
import sys
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

def create_presentation():
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    # Color Palette (Dark Theme)
    BG_COLOR = RGBColor(15, 23, 42)       # #0f172a
    CARD_BG = RGBColor(30, 41, 59)        # #1e293b
    ACCENT_BLUE = RGBColor(56, 189, 248)  # #38bdf8
    ACCENT_INDIGO = RGBColor(129, 140, 248)# #818cf8
    TEXT_MAIN = RGBColor(248, 250, 252)   # #f8fafc
    TEXT_MUTED = RGBColor(148, 163, 184)  # #94a3b8

    blank_layout = prs.slide_layouts[6]

    def add_slide_background(slide):
        bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, prs.slide_width, prs.slide_height)
        bg.fill.solid()
        bg.fill.fore_color.rgb = BG_COLOR
        bg.line.fill.background()
        return bg

    def add_header(slide, title_text, subtitle_text="THz-Unified-Optimizer | Еженедельный отчёт"):
        txBox = slide.shapes.add_textbox(Inches(0.8), Inches(0.4), Inches(11.7), Inches(1.1))
        tf = txBox.text_frame
        tf.word_wrap = True
        
        p = tf.paragraphs[0]
        p.text = title_text
        p.font.size = Pt(28)
        p.font.bold = True
        p.font.color.rgb = ACCENT_BLUE
        p.font.name = 'Segoe UI'

        if subtitle_text:
            p2 = tf.add_paragraph()
            p2.text = subtitle_text
            p2.font.size = Pt(14)
            p2.font.color.rgb = TEXT_MUTED
            p2.font.name = 'Segoe UI'

    # ------------------ SLIDE 1: Title ------------------
    slide1 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide1)

    tb = slide1.shapes.add_textbox(Inches(1.0), Inches(2.0), Inches(11.3), Inches(3.5))
    tf = tb.text_frame
    tf.word_wrap = True

    p = tf.paragraphs[0]
    p.text = "THz-Unified-Optimizer"
    p.font.size = Pt(44)
    p.font.bold = True
    p.font.color.rgb = ACCENT_BLUE

    p2 = tf.add_paragraph()
    p2.text = "Результаты недели (14 – 21 июля 2026)"
    p2.font.size = Pt(28)
    p2.font.color.rgb = ACCENT_INDIGO

    p3 = tf.add_paragraph()
    p3.text = "\nПрецизионное физико-статистическое моделирование проволочных решеток (WGP),\nанализ невязок, мультистартовая оптимизация и валидация на образце 40/20 мкм."
    p3.font.size = Pt(16)
    p3.font.color.rgb = TEXT_MAIN

    p4 = tf.add_paragraph()
    p4.text = "\nАвтор: Antigravity IDE | Лаборатория ТГц спектроскопии"
    p4.font.size = Pt(14)
    p4.font.color.rgb = TEXT_MUTED

    # ------------------ SLIDE 2: Key Milestones ------------------
    slide2 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide2)
    add_header(slide2, "🎯 Главные вехи и достижения недели")

    milestones = [
        ("Физическая эволюция модели", "Переход от Pure Blanco к комплексной Jones-матрице с учетом поверхностного импеданса Друде (вольфрам), рассеяния e^(-α ν^γ) и фазовой задержки τ_par."),
        ("Статистический анализ & LMFIT", "Построение матриц ковариации, доверительных интервалов (C(P,D) = -0.9991) и статистический доказательный анализ через Ablation Study (AIC/BIC)."),
        ("Анализ невязок (AGENTS Rule #5)", "Подтверждение ненормальности остатков (Shapiro-Wilk) и 85% кросс-серийной корреляции при скользящих углах 50°–90° (апертурная дифракция)."),
        ("Эксперимент на решетке 40/20 мкм", "Измерение серий grid_40_20, подтверждение закона масштабирования Deff = Dphys * (1 - 0.85 * D/P) и снижение χ²_ν до 0.00310.")
    ]

    for i, (m_title, m_desc) in enumerate(milestones):
        top_pos = Inches(1.6 + i * 1.3)
        shape = slide2.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(0.8), top_pos, Inches(11.733), Inches(1.1))
        shape.fill.solid()
        shape.fill.fore_color.rgb = CARD_BG
        shape.line.color.rgb = ACCENT_INDIGO

        tf = shape.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        p.text = f"{i+1}. {m_title}"
        p.font.size = Pt(18)
        p.font.bold = True
        p.font.color.rgb = ACCENT_BLUE

        p2 = tf.add_paragraph()
        p2.text = m_desc
        p2.font.size = Pt(13)
        p2.font.color.rgb = TEXT_MAIN

    # ------------------ SLIDE 3: Ablation Study ------------------
    slide3 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide3)
    add_header(slide3, "📈 Модельный Ablation Study (Сравнение AIC/BIC)")

    # Left: Table
    rows, cols = 5, 6
    left, top, width, height = Inches(0.8), Inches(1.8), Inches(6.8), Inches(4.5)
    table_shape = slide3.shapes.add_table(rows, cols, left, top, width, height)
    table = table_shape.table

    headers = ["Модель", "D_eff (мкм)", "χ²_ν", "RMSE (дБ)", "AIC", "ΔAIC"]
    data = [
        ["M0: Blanco", "4.614", "0.00149", "0.039", "−19655.7", "0.0"],
        ["M1: +Drude", "4.462", "0.00178", "0.042", "−19116.0", "+539.8 ❌"],
        ["M2: +Scat (γ=2)", "4.394", "0.00157", "0.040", "−19491.2", "−375.2 ✅"],
        ["M3: +Scat (free γ)", "4.384", "0.00157", "0.040", "−19492.4", "−1.2 ≈0"]
    ]

    for col_idx, h in enumerate(headers):
        cell = table.cell(0, col_idx)
        cell.fill.solid()
        cell.fill.fore_color.rgb = CARD_BG
        p = cell.text_frame.paragraphs[0]
        p.text = h
        p.font.bold = True
        p.font.size = Pt(12)
        p.font.color.rgb = ACCENT_BLUE

    for row_idx, r_data in enumerate(data):
        for col_idx, val in enumerate(r_data):
            cell = table.cell(row_idx + 1, col_idx)
            cell.fill.solid()
            cell.fill.fore_color.rgb = CARD_BG
            p = cell.text_frame.paragraphs[0]
            p.text = val
            p.font.size = Pt(11)
            p.font.color.rgb = TEXT_MAIN

    # Right: Image if available
    img_path = r"c:\THz-Unified-Optimizer\docs\images\model_ablation_comparison.png"
    if os.path.exists(img_path):
        slide3.shapes.add_picture(img_path, Inches(7.8), Inches(1.8), Inches(4.7), Inches(4.5))

    # ------------------ SLIDE 4: 40/20 Sample & Scaling Law ------------------
    slide4 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide4)
    add_header(slide4, "🧪 Новая решётка 40/20 мкм и закон масштабирования")

    img_scaling = r"c:\THz-Unified-Optimizer\docs\images\deff_scaling_law.png"
    if os.path.exists(img_scaling):
        slide4.shapes.add_picture(img_scaling, Inches(0.8), Inches(1.8), Inches(5.5), Inches(4.8))

    # Text box on right
    tb = slide4.shapes.add_textbox(Inches(6.6), Inches(1.8), Inches(5.9), Inches(4.8))
    tf = tb.text_frame
    tf.word_wrap = True

    p = tf.paragraphs[0]
    p.text = "Экспериментальная проверка геометрического блокирования"
    p.font.size = Pt(20)
    p.font.bold = True
    p.font.color.rgb = ACCENT_BLUE

    points = [
        "Старая решётка (P=15.5 мкм, D=5.67 мкм, D/P=0.71): Deff/Dphys = 0.40 (Deff = 4.38 мкм).",
        "Новая решётка (P=40.0 мкм, D=20.0 мкм, D/P=0.50): Deff/Dphys = 0.57 (Deff = 11.38 мкм).",
        "Сформулирован эмпирический закон: Deff = Dphys * (1 - 0.85 * D/P).",
        "Подтверждение: при уменьшении плотности WGP effective diameter стремится к физическому полосовому диаметру Dphys."
    ]
    for pt in points:
        p_pt = tf.add_paragraph()
        p_pt.text = f"• {pt}"
        p_pt.font.size = Pt(14)
        p_pt.font.color.rgb = TEXT_MAIN

    # ------------------ SLIDE 5: Residuals & Phase Delay ------------------
    slide5 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide5)
    add_header(slide5, "🔬 Анализ невязок и фазовая анизотропия τ_par")

    img_res = r"c:\THz-Unified-Optimizer\docs\images\residuals_comprehensive_maps.png"
    if os.path.exists(img_res):
        slide5.shapes.add_picture(img_res, Inches(0.8), Inches(1.8), Inches(5.6), Inches(4.8))

    img_td = r"c:\THz-Unified-Optimizer\docs\images\time_delta_vs_angle_40_20.png"
    if os.path.exists(img_td):
        slide5.shapes.add_picture(img_td, Inches(6.7), Inches(1.8), Inches(5.8), Inches(4.8))

    # ------------------ SLIDE 6: Conclusions & Next Steps ------------------
    slide6 = prs.slides.add_slide(blank_layout)
    add_slide_background(slide6)
    add_header(slide6, "🎯 Выводы и дальнейшие шаги")

    tb = slide6.shapes.add_textbox(Inches(0.8), Inches(1.8), Inches(5.6), Inches(4.8))
    tf = tb.text_frame
    tf.word_wrap = True

    p = tf.paragraphs[0]
    p.text = "✅ Итоги анализа текущих данных"
    p.font.size = Pt(20)
    p.font.bold = True
    p.font.color.rgb = ACCENT_BLUE

    items_left = [
        "Аналитический потенциал модели Бланко-Друде исчерпан полностью (RMSE = 0.427 дБ).",
        "Точность вычисления параметров достигла <0.2% (Deff = 4.820 ± 0.007 мкм).",
        "Дифракция на оправе ложна (апертура 4\" = 101.6 мм при пучке 12 мм). Невязка 85% на 50°–90° связана с ограничением С/Ш детектора при глубоком затухании.",
        "Исключены сбойные серии с асимметрией >70%."
    ]
    for it in items_left:
        p_it = tf.add_paragraph()
        p_it.text = f"• {it}"
        p_it.font.size = Pt(14)
        p_it.font.color.rgb = TEXT_MAIN

    tb2 = slide6.shapes.add_textbox(Inches(6.8), Inches(1.8), Inches(5.7), Inches(4.8))
    tf2 = tb2.text_frame
    tf2.word_wrap = True

    p_r = tf2.paragraphs[0]
    p_r.text = "🚀 План будущих работ"
    p_r.font.size = Pt(20)
    p_r.font.bold = True
    p_r.font.color.rgb = ACCENT_INDIGO

    items_right = [
        "FTIR Спектрометрия (Bruker Vertex 70): съем данных в диапазоне > 3 ТГц.",
        "3D методы (RCWA / FDTD): переподготовка модели для голографических пленок (1200 лин/мм).",
        "Журнальная публикация: подготовка статьи по соотношениям Deff(D/P) и анизотропной фазовой задержке."
    ]
    for it in items_right:
        p_it = tf2.add_paragraph()
        p_it.text = f"• {it}"
        p_it.font.size = Pt(14)
        p_it.font.color.rgb = TEXT_MAIN

    out_path = r"c:\THz-Unified-Optimizer\docs\artifacts\presentation_weekly.pptx"
    prs.save(out_path)
    print(f"Presentation saved successfully to {out_path}")

if __name__ == '__main__':
    create_presentation()
