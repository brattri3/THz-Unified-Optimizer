import os
import sys
from pathlib import Path
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

def main():
    pdf_path = Path("c:/THz-Unified-Optimizer/docs/artifacts/final_comprehensive_report.pdf")
    pdf_path.parent.mkdir(exist_ok=True, parents=True)
    
    # Register Cyrillic font (Arial is standard on Windows)
    arial_path = "C:/Windows/Fonts/arial.ttf"
    arial_bd_path = "C:/Windows/Fonts/arialbd.ttf"
    
    if os.path.exists(arial_path):
        pdfmetrics.registerFont(TTFont('Arial', arial_path))
    else:
        # Fallback to Helvetica if Arial is not found (might not render Cyrillic correctly, but avoids crash)
        print("Arial font not found, Cyrillic text may not render correctly.")
        
    if os.path.exists(arial_bd_path):
        pdfmetrics.registerFont(TTFont('Arial-Bold', arial_bd_path))
        
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        rightMargin=40, leftMargin=40, topMargin=40, bottomMargin=40
    )
    
    styles = getSampleStyleSheet()
    
    # Create custom styles supporting Russian language (Arial font)
    title_style = ParagraphStyle(
        'TitleStyle',
        parent=styles['Heading1'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=20,
        leading=24,
        textColor=colors.HexColor('#1A365D'),
        spaceAfter=15,
        alignment=1 # Center
    )
    
    subtitle_style = ParagraphStyle(
        'SubtitleStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=11,
        leading=14,
        textColor=colors.HexColor('#4A5568'),
        spaceAfter=25,
        alignment=1 # Center
    )
    
    h1_style = ParagraphStyle(
        'H1Style',
        parent=styles['Heading2'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=14,
        leading=18,
        textColor=colors.HexColor('#2B6CB0'),
        spaceBefore=12,
        spaceAfter=8
    )
    
    h2_style = ParagraphStyle(
        'H2Style',
        parent=styles['Heading3'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=12,
        leading=15,
        textColor=colors.HexColor('#2D3748'),
        spaceBefore=8,
        spaceAfter=6
    )
    
    body_style = ParagraphStyle(
        'BodyStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=10,
        leading=14,
        textColor=colors.HexColor('#2D3748'),
        spaceAfter=10
    )
    
    table_text_style = ParagraphStyle(
        'TableText',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=9,
        leading=11
    )
    
    table_header_style = ParagraphStyle(
        'TableHeader',
        parent=styles['Normal'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=9,
        leading=11
    )

    story = []
    
    # Title Page / Header
    story.append(Paragraph("Полный аналитический отчет:<br/>Моделирование и анализ невязок THz-TDS", title_style))
    story.append(Paragraph("Исследование угловых зависимостей, комплексная оценка невязок и кросс-корреляционный анализ экспериментальных данных", subtitle_style))
    story.append(Spacer(1, 10))
    
    # Part 1
    story.append(Paragraph("Часть 1. Описание физической модели Бланко и параметры", h1_style))
    story.append(Paragraph(
        "В данном исследовании представлено сопоставление экспериментальных данных ТГц-TDS спектроскопии с теоретической электродинамической моделью Бланко для системы из двух проволочных поляризаторов (аттенюатора). "
        "В отличие от базовой версии модели, мы внедрили <b>физическую модель импеданса Друде</b> (для учета частотно-зависимых омических потерь в вольфраме) и <b>фазовую задержку оптического пути</b>.",
        body_style
    ))
    story.append(Paragraph("<b>Оптимизированные параметры решётки и тракта:</b>", body_style))
    story.append(Paragraph("- Период проволочной решётки ($P$): 15.500 мкм (зафиксирован на паспортном значении)", body_style))
    story.append(Paragraph("- Диаметр проволоки ($D$): 4.398 мкм", body_style))
    story.append(Paragraph("- Систематический сдвиг угла (theta_offset): -0.05 deg", body_style))
    story.append(Paragraph("- Коэффициент рассеяния (loss_factor): 0.316 (со степенью gamma = 1.06)", body_style))
    story.append(Paragraph("- Фазовая задержка (tau_ps): 0.033 пс", body_style))
    story.append(Spacer(1, 10))
    
    # Part 2
    story.append(Paragraph("Часть 2. Прямое сравнение: Спектральный и Интегральный методы", h1_style))
    story.append(Paragraph(
        "Амплитудный коэффициент пропускания оценивался двумя способами:<br/>"
        "1. <b>Спектральный анализ:</b> T(nu) = |Es(nu)| / |Eb(nu)| на фиксированных частотах 0.5 ТГц и 1.0 ТГц.<br/>"
        "2. <b>Интегральный метод:</b> Оценка по полной энергии импульса во временной области (с учетом RMS-усреднения теоретической модели по фоновому спектру).",
        body_style
    ))
    
    story.append(Paragraph("<b>Серия измерений 356att</b>", h2_style))
    img_spec_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_356att.png")
    img_int_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_356att.png")
    
    if img_spec_356.exists():
        story.append(Image(str(img_spec_356), width=450, height=220))
        story.append(Spacer(1, 5))
    if img_int_356.exists():
        story.append(Image(str(img_int_356), width=450, height=220))
        story.append(Spacer(1, 10))
        
    story.append(PageBreak())
    
    story.append(Paragraph("<b>Серия измерений series3</b>", h2_style))
    img_spec_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_series3.png")
    img_int_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_series3.png")
    
    if img_spec_s3.exists():
        story.append(Image(str(img_spec_s3), width=450, height=220))
        story.append(Spacer(1, 5))
    if img_int_s3.exists():
        story.append(Image(str(img_int_s3), width=450, height=220))
        story.append(Spacer(1, 10))
        
    # Part 3
    story.append(Paragraph("Часть 3. Комплексный анализ невязок (Амплитуда и Фаза)", h1_style))
    story.append(Paragraph(
        "Для оценки качества подгонки комплексное отношение Es(nu)/Eb(nu) было разделено на амплитудную (delta A, дБ) и фазовую (delta phi, рад) невязки с исключением линий поглощения водяного пара.",
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
    t = Table(table_data, colWidths=[200, 150, 150])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E2E8F0')),
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('BOTTOMPADDING', (0,0), (-1,0), 6),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
    ]))
    story.append(t)
    story.append(Spacer(1, 15))
    
    img_res_amp_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_amp_356att.png")
    img_res_phase_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_phase_356att.png")
    
    if img_res_amp_356.exists():
        story.append(Paragraph("<b>2D Карты невязок для 356att (Амплитуда и Фаза)</b>", h2_style))
        story.append(Image(str(img_res_amp_356), width=450, height=180))
        story.append(Spacer(1, 5))
    if img_res_phase_356.exists():
        story.append(Image(str(img_res_phase_356), width=450, height=180))
        story.append(Spacer(1, 10))
        
    story.append(PageBreak())
    
    img_res_amp_s3 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_amp_series3.png")
    img_res_phase_s3 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_phase_series3.png")
    
    if img_res_amp_s3.exists():
        story.append(Paragraph("<b>2D Карты невязок для series3 (Амплитуда и Фаза)</b>", h2_style))
        story.append(Image(str(img_res_amp_s3), width=450, height=180))
        story.append(Spacer(1, 5))
    if img_res_phase_s3.exists():
        story.append(Image(str(img_res_phase_s3), width=450, height=180))
        story.append(Spacer(1, 10))
        
    # Part 4
    story.append(Paragraph("Часть 4. Корреляционный анализ отклонений (356att vs series3)", h1_style))
    story.append(Paragraph(
        "Для проверки гипотезы о систематической природе оставшихся отклонений мы провели кросс-корреляционный анализ невязок между независимыми сериями 356att и series3 на 12 совпадающих угловых позициях (1357 общих частотно-угловых точек).",
        body_style
    ))
    story.append(Paragraph("<b>Глобальная корреляция:</b>", body_style))
    story.append(Paragraph("- Амплитуда: R_A = 0.282 (p-value < 0.001)", body_style))
    story.append(Paragraph("- Фаза: R_phi = 0.594 (p-value < 0.001)", body_style))
    
    img_corr_ang = Path("c:/THz-Unified-Optimizer/docs/images/correlation_angular.png")
    img_corr_scat = Path("c:/THz-Unified-Optimizer/docs/images/correlation_scatter.png")
    
    if img_corr_ang.exists():
        story.append(Image(str(img_corr_ang), width=450, height=220))
        story.append(Spacer(1, 5))
    if img_corr_scat.exists():
        story.append(Image(str(img_corr_scat), width=450, height=200))
        story.append(Spacer(1, 10))
        
    # Part 5
    story.append(Paragraph("Часть 5. Итоговое физическое заключение", h1_style))
    story.append(Paragraph(
        "Модель Blanco + Drude достигла своего теоретического предела для данного набора данных (ошибка < 1 дБ). "
        "Оставшиеся невязки — это не случайный шум измерений, а реальные физические или геометрические систематические эффекты экспериментального стенда. К ним относятся:<br/>"
        "- <b>Дифракция на краях апертуры</b> (пучок перестает быть параксиальным и задевает оправку при больших углах вращения).<br/>"
        "- <b>Оптическая анизотропия</b> подложки или локальная неидеальность натяжения нитей вольфрамовой сетки.",
        body_style
    ))
    
    doc.build(story)
    print("PDF successfully generated.")

if __name__ == "__main__":
    main()
