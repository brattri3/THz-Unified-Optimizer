import os
import sys
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

def make_side_by_side_images(img_path1, img_path2, width=240, height=120):
    if img_path1.exists() and img_path2.exists():
        im1 = Image(str(img_path1), width=width, height=height)
        im2 = Image(str(img_path2), width=width, height=height)
        t = Table([[im1, im2]], colWidths=[width + 10, width + 10])
        t.setStyle(TableStyle([
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
            ('BOTTOMPADDING', (0,0), (-1,-1), 0),
            ('TOPPADDING', (0,0), (-1,-1), 0),
            ('LEFTPADDING', (0,0), (-1,-1), 0),
            ('RIGHTPADDING', (0,0), (-1,-1), 0),
        ]))
        return t
    elif img_path1.exists():
        return Image(str(img_path1), width=width*1.5, height=height*1.5)
    elif img_path2.exists():
        return Image(str(img_path2), width=width*1.5, height=height*1.5)
    return None

def main():
    pdf_path = Path("c:/THz-Unified-Optimizer/docs/artifacts/final_comprehensive_report.pdf")
    pdf_path.parent.mkdir(exist_ok=True, parents=True)
    
    # Register Cyrillic font
    arial_path = "C:/Windows/Fonts/arial.ttf"
    arial_bd_path = "C:/Windows/Fonts/arialbd.ttf"
    
    if os.path.exists(arial_path):
        pdfmetrics.registerFont(TTFont('Arial', arial_path))
    else:
        print("Arial font not found, Cyrillic text may not render correctly.")
        
    if os.path.exists(arial_bd_path):
        pdfmetrics.registerFont(TTFont('Arial-Bold', arial_bd_path))
        
    # Standard academic paper margins: 0.75 in (54 pt) or 0.5 in (36 pt)
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        rightMargin=36, leftMargin=36, topMargin=40, bottomMargin=40
    )
    
    styles = getSampleStyleSheet()
    
    # Academic journal styling: Times/Arial, single-spaced, smaller headings
    title_style = ParagraphStyle(
        'TitleStyle',
        parent=styles['Heading1'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=16,
        leading=20,
        textColor=colors.HexColor('#1A202C'),
        spaceAfter=10,
        alignment=1
    )
    
    subtitle_style = ParagraphStyle(
        'SubtitleStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=10,
        leading=13,
        textColor=colors.HexColor('#4A5568'),
        spaceAfter=15,
        alignment=1
    )
    
    h1_style = ParagraphStyle(
        'H1Style',
        parent=styles['Heading2'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=12,
        leading=15,
        textColor=colors.HexColor('#2D3748'),
        spaceBefore=12,
        spaceAfter=6,
        keepWithNext=True
    )
    
    h2_style = ParagraphStyle(
        'H2Style',
        parent=styles['Heading3'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=10,
        leading=13,
        textColor=colors.HexColor('#2D3748'),
        spaceBefore=8,
        spaceAfter=4,
        keepWithNext=True
    )
    
    body_style = ParagraphStyle(
        'BodyStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=9.5,
        leading=13,
        textColor=colors.HexColor('#2D3748'),
        spaceAfter=8,
        alignment=4 # Justified
    )
    
    table_text_style = ParagraphStyle(
        'TableText',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=8.5,
        leading=11
    )
    
    table_header_style = ParagraphStyle(
        'TableHeader',
        parent=styles['Normal'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=8.5,
        leading=11
    )
    
    caption_style = ParagraphStyle(
        'CaptionStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=8,
        leading=10,
        textColor=colors.HexColor('#4A5568'),
        alignment=1,
        spaceAfter=10
    )

    story = []
    
    # Document Header
    story.append(Paragraph("Полный аналитический отчет:<br/>Моделирование и анализ невязок THz-TDS", title_style))
    story.append(Paragraph("Исследование угловых зависимостей, комплексная оценка невязок и кросс-корреляционный анализ экспериментальных данных", subtitle_style))
    story.append(Spacer(1, 5))
    
    # Part 1
    story.append(Paragraph("1. Описание физической модели Бланко и параметры", h1_style))
    story.append(Paragraph(
        "В данном исследовании представлено сопоставление экспериментальных данных ТГц-TDS спектроскопии с теоретической электродинамической моделью Бланко для системы из двух проволочных поляризаторов (аттенюатора). "
        "В отличие от базовой версии модели, мы внедрили <b>физическую модель импеданса Друде</b> (для учета частотно-зависимых омических потерь в вольфраме) и <b>фазовую задержку оптического пути</b>.",
        body_style
    ))
    story.append(Paragraph(
        "<b>Оптимизированные параметры решётки и тракта (Global Average):</b><br/>"
        "• Период проволочной решётки (<i>P</i>): 15.500 мкм (зафиксирован на паспортном значении)<br/>"
        "• Диаметр проволоки (<i>D</i>): 4.398 мкм<br/>"
        "• Систематический сдвиг угла (&theta;<sub>offset</sub>): -0.05°<br/>"
        "• Коэффициент рассеяния (<i>loss_factor</i>): 0.316 (со степенью &gamma; = 1.06)<br/>"
        "• Фазовая задержка (&tau;<sub>ps</sub>): 0.033 пс",
        body_style
    ))
    
    # Part 2
    story.append(Paragraph("2. Прямое сравнение: Спектральный и Интегральный методы", h1_style))
    story.append(Paragraph(
        "Амплитудный коэффициент пропускания оценивался двумя способами:<br/>"
        "1. <b>Спектральный анализ:</b> T(&nu;) = |Es(&nu;)| / |Eb(&nu;)| на фиксированных частотах 0.5 ТГц и 1.0 ТГц.<br/>"
        "2. <b>Интегральный метод:</b> Оценка по полной энергии импульса во временной области (с учетом RMS-усреднения теоретической модели по фоновому спектру).",
        body_style
    ))
    
    # Images for 356att side-by-side
    story.append(Paragraph("<b>Серия измерений 356att</b>", h2_style))
    img_spec_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_356att.png")
    img_int_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_356att.png")
    fig_356 = make_side_by_side_images(img_spec_356, img_int_356, width=250, height=125)
    if fig_356:
        story.append(fig_356)
        story.append(Paragraph("Рис 1. Сравнение экспериментальных точек с теорией для 356att: спектральный анализ (слева) и интегральный метод (справа).", caption_style))
        
    # Images for series3 side-by-side
    story.append(Paragraph("<b>Серия измерений series3</b>", h2_style))
    img_spec_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_series3.png")
    img_int_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_series3.png")
    fig_s3 = make_side_by_side_images(img_spec_s3, img_int_s3, width=250, height=125)
    if fig_s3:
        story.append(fig_s3)
        story.append(Paragraph("Рис 2. Сравнение экспериментальных точек с теорией для series3: спектральный анализ (слева) и интегральный метод (справа).", caption_style))
        
    # Part 3
    story.append(Paragraph("3. Комплексный анализ невязок (Амплитуда и Фаза)", h1_style))
    story.append(Paragraph(
        "Для оценки качества подгонки комплексное отношение Es(&nu;)/Eb(&nu;) было разделено на амплитудную (&Delta;A, дБ) и фазовую (&Delta;&phi;, рад) невязки с исключением линий поглощения водяного пара.",
        body_style
    ))
    
    # Table of stats
    table_data = [
        [Paragraph("<b>Серия / Метрика</b>", table_header_style), Paragraph("<b>Шапиро-Уилк (p-value)</b>", table_header_style), Paragraph("<b>Харке-Бера (p-value)</b>", table_header_style)],
        [Paragraph("356att (Амплитуда)", table_text_style), Paragraph("8.42e-36", table_text_style), Paragraph("4.51e-212", table_text_style)],
        [Paragraph("356att (Фаза)", table_text_style), Paragraph("2.35e-39", table_text_style), Paragraph("9.13e-261", table_text_style)],
        [Paragraph("series3 (Амплитуда)", table_text_style), Paragraph("8.44e-17", table_text_style), Paragraph("6.15e-27", table_text_style)],
        [Paragraph("series3 (Фаза)", table_text_style), Paragraph("1.75e-52", table_text_style), Paragraph("0.00", table_text_style)]
    ]
    t = Table(table_data, colWidths=[180, 160, 160])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#EDF2F7')),
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
    ]))
    story.append(t)
    story.append(Spacer(1, 8))
    
    # Residuals images side-by-side for 356att
    story.append(Paragraph("<b>2D Карты невязок для 356att (Амплитуда и Фаза)</b>", h2_style))
    img_res_amp_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_amp_356att.png")
    img_res_phase_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_phase_356att.png")
    fig_res_356 = make_side_by_side_images(img_res_amp_356, img_res_phase_356, width=250, height=110)
    if fig_res_356:
        story.append(fig_res_356)
        story.append(Paragraph("Рис 3. Распределение амплитудных (слева) и фазовых (справа) невязок по частотам и углам для серии 356att.", caption_style))
        
    # Residuals images side-by-side for series3
    story.append(Paragraph("<b>2D Карты невязок для series3 (Амплитуда и Фаза)</b>", h2_style))
    img_res_amp_s3 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_amp_series3.png")
    img_res_phase_s3 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_phase_series3.png")
    fig_res_s3 = make_side_by_side_images(img_res_amp_s3, img_res_phase_s3, width=250, height=110)
    if fig_res_s3:
        story.append(fig_res_s3)
        story.append(Paragraph("Рис 4. Распределение амплитудных (слева) и фазовых (справа) невязок по частотам и углам для серии series3.", caption_style))
        
    # Part 4
    story.append(Paragraph("4. Корреляционный анализ отклонений (356att vs series3)", h1_style))
    story.append(Paragraph(
        "Для проверки гипотезы о систематической природе оставшихся отклонений мы провели кросс-корреляционный анализ невязок между независимыми сериями 356att и series3 на 12 совпадающих угловых позициях (1357 общих частотно-угловых точек).",
        body_style
    ))
    story.append(Paragraph(
        "<b>Глобальная корреляция:</b><br/>"
        "• Амплитуда: R_A = 0.282 (p-value < 0.001)<br/>"
        "• Фаза: R_phi = 0.594 (p-value < 0.001)",
        body_style
    ))
    
    img_corr_ang = Path("c:/THz-Unified-Optimizer/docs/images/correlation_angular.png")
    img_corr_scat = Path("c:/THz-Unified-Optimizer/docs/images/correlation_scatter.png")
    fig_corr = make_side_by_side_images(img_corr_ang, img_corr_scat, width=250, height=120)
    if fig_corr:
        story.append(fig_corr)
        story.append(Paragraph("Рис 5. Угловая зависимость корреляции невязок (слева) и глобальная корреляция амплитуды/фазы (справа).", caption_style))
        
    # Part 5
    story.append(Paragraph("5. Итоговое физическое заключение", h1_style))
    story.append(Paragraph(
        "Модель Blanco + Drude достигла своего теоретического предела для данного набора данных (ошибка < 1 дБ). "
        "Оставшиеся невязки — это не случайный шум измерений, а реальные физические или геометрические систематические эффекты экспериментального стенда. К ним относятся:<br/>"
        "• <b>Дифракция на краях апертуры</b> (пучок перестает быть параксиальным и задевает оправку при больших углах вращения).<br/>"
        "• <b>Оптическая анизотропия</b> подложки или локальная неидеальность натяжения нитей вольфрамовой сетки.",
        body_style
    ))
    
    doc.build(story)
    print("Academic PDF successfully generated.")

if __name__ == "__main__":
    main()
