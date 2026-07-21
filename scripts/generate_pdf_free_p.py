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
    pdf_path = Path("c:/THz-Unified-Optimizer/docs/artifacts/final_comprehensive_report_free_p.pdf")
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
        
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        rightMargin=36, leftMargin=36, topMargin=40, bottomMargin=40
    )
    
    styles = getSampleStyleSheet()
    
    title_style = ParagraphStyle(
        'TitleStyle',
        parent=styles['Heading1'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=15,
        leading=18,
        textColor=colors.HexColor('#1A202C'),
        spaceAfter=8,
        alignment=1
    )
    
    subtitle_style = ParagraphStyle(
        'SubtitleStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=9.5,
        leading=12,
        textColor=colors.HexColor('#4A5568'),
        spaceAfter=12,
        alignment=1
    )
    
    h1_style = ParagraphStyle(
        'H1Style',
        parent=styles['Heading2'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=11.5,
        leading=14,
        textColor=colors.HexColor('#2D3748'),
        spaceBefore=10,
        spaceAfter=5,
        keepWithNext=True
    )
    
    h2_style = ParagraphStyle(
        'H2Style',
        parent=styles['Heading3'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=9.5,
        leading=12,
        textColor=colors.HexColor('#2D3748'),
        spaceBefore=6,
        spaceAfter=3,
        keepWithNext=True
    )
    
    body_style = ParagraphStyle(
        'BodyStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=9,
        leading=12,
        textColor=colors.HexColor('#2D3748'),
        spaceAfter=6,
        alignment=4 # Justified
    )
    
    table_text_style = ParagraphStyle(
        'TableText',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=8,
        leading=10
    )
    
    table_header_style = ParagraphStyle(
        'TableHeader',
        parent=styles['Normal'],
        fontName='Arial-Bold' if os.path.exists(arial_bd_path) else 'Arial',
        fontSize=8,
        leading=10
    )
    
    caption_style = ParagraphStyle(
        'CaptionStyle',
        parent=styles['Normal'],
        fontName='Arial',
        fontSize=7.5,
        leading=9,
        textColor=colors.HexColor('#4A5568'),
        alignment=1,
        spaceAfter=8
    )

    story = []
    
    # Document Header
    story.append(Paragraph("Аналитический отчет (Свободный период): Моделирование и анализ погрешностей THz-TDS", title_style))
    story.append(Paragraph("Результаты комплексной оптимизации модели Бланко-Друде со свободным шагом P и оценкой ошибок методом Гессиана", subtitle_style))
    story.append(Spacer(1, 4))
    
    # Part 1
    story.append(Paragraph("1. Введение и постановка задачи", h1_style))
    story.append(Paragraph(
        "В данном исследовании представлен анализ комплексной подгонки Blanco-Drude модели со свободным периодом решетки $P$ и оценка неопределенностей параметров с использованием численного Гессиана. "
        "Период оптимизировался в границах [14.5, 17.5] мкм, а диаметр — [4.8, 6.2] мкм. Оценка погрешностей $\pm \sigma$ проводилась на основе диагональных элементов ковариационной матрицы $\Sigma = \\frac{2 L_{min}}{M} H^{-1}$ в точке минимума.",
        body_style
    ))
    
    # Part 2
    story.append(Paragraph("2. Результаты оптимизации и оценки Гессиана", h1_style))
    
    # Parameters Table
    table_data = [
        [Paragraph("<b>Параметр</b>", table_header_style), Paragraph("<b>Серия 356att (M=2583)</b>", table_header_style), Paragraph("<b>Серия series3 (M=1388)</b>", table_header_style), Paragraph("<b>Global Average (M=2000)</b>", table_header_style)],
        [Paragraph("Период P (мкм)", table_text_style), Paragraph("16.729 ± 0.025", table_text_style), Paragraph("17.500 ± 0.000 (границ.)", table_text_style), Paragraph("17.500 ± 0.000 (границ.)", table_text_style)],
        [Paragraph("Диаметр D (мкм)", table_text_style), Paragraph("4.820 ± 0.007", table_text_style), Paragraph("4.991 ± 0.020", table_text_style), Paragraph("5.053 ± 0.010", table_text_style)],
        [Paragraph("Шумовой factor", table_text_style), Paragraph("0.287 ± 0.0016", table_text_style), Paragraph("0.0017 ± 0.0245", table_text_style), Paragraph("0.287 ± 0.0023", table_text_style)],
        [Paragraph("Показатель степени gamma", table_text_style), Paragraph("0.294 ± 0.0016", table_text_style), Paragraph("0.721 ± 0.000", table_text_style), Paragraph("1.868 ± 0.0042", table_text_style)],
        [Paragraph("Сдвиг угла theta_offset (deg)", table_text_style), Paragraph("-0.013 ± 0.0016", table_text_style), Paragraph("0.004 ± 0.0825", table_text_style), Paragraph("-0.183 ± 0.0023", table_text_style)],
        [Paragraph("Задержка фазы tau_ps (ps)", table_text_style), Paragraph("0.0237 ± 0.0008", table_text_style), Paragraph("0.0125 ± 0.0009", table_text_style), Paragraph("0.0327 ± 0.0010", table_text_style)],
        [Paragraph("Минимум Loss (L_min)", table_text_style), Paragraph("0.00273", table_text_style), Paragraph("0.00194", table_text_style), Paragraph("0.00266", table_text_style)]
    ]
    
    t = Table(table_data, colWidths=[140, 120, 120, 120])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#EDF2F7')),
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3),
        ('TOPPADDING', (0,0), (-1,-1), 3),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E0')),
    ]))
    story.append(t)
    story.append(Spacer(1, 8))
    
    # Part 3
    story.append(Paragraph("3. Физическое обсуждение погрешностей и чувствительности", h1_style))
    story.append(Paragraph(
        "<b>3.1. Высокая геометрическая чувствительность:</b><br/>"
        "Для серии 356att погрешности определения геометрических параметров крайне малы: &sigma;_P = ±25 нм и &sigma;_D = ±7 нм. "
        "Это доказывает, что целевая функция имеет четко локализованный глобальный минимум в P-D пространстве. "
        "Модель Blanco-Drude чрезвычайно чувствительна к изменению пространственной геометрии на уровне десятков нанометров, что подтверждает физическую адекватность ее фазового и амплитудного отклика.",
        body_style
    ))
    story.append(Paragraph(
        "<b>3.2. Выход периода P на верхнюю границу (17.5 мкм):</b><br/>"
        "Для серии series3 и Global Average период P сошелся ровно к верхнему лимиту 17.500 мкм (с нулевой погрешностью из-за вырождения штрафной функции). "
        "Это указывает на физическое различие между амплитудными и фазовыми спектрами: для согласования малого поглощения (loss_factor ~ 0.0017 в series3) и фазовой задержки модель вынуждена увеличивать шаг решетки. "
        "Вероятная причина — <i>непараксиальность сфокусированного пучка (наклонное падение)</i>, из-за чего проекционный период P' = P/cos(&theta;_inc) превышает паспортный.",
        body_style
    ))
    story.append(Paragraph(
        "<b>3.3. Точность фазового сдвига:</b><br/>"
        "Погрешность фазовой задержки &tau;_ps составила всего ±0.8--1.0 фс (фемтосекунды). "
        "Это демонстрирует высокую метрологическую точность ТГц-TDS приборов при фазовом фитинге и доказывает важность включения фазовой задержки в модель.",
        body_style
    ))
    
    # Part 4
    story.append(Paragraph("4. Прямое сравнение и невязки", h1_style))
    
    # 356att side-by-side
    story.append(Paragraph("<b>Серия измерений 356att</b>", h2_style))
    img_spec_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_356att.png")
    img_int_356 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_356att.png")
    fig_356 = make_side_by_side_images(img_spec_356, img_int_356, width=245, height=120)
    if fig_356:
        story.append(fig_356)
        story.append(Paragraph("Рис 1. 356att: спектральный анализ пропускания (слева) и интегральный метод (справа).", caption_style))
        
    # series3 side-by-side
    story.append(Paragraph("<b>Серия измерений series3</b>", h2_style))
    img_spec_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_spectral_series3.png")
    img_int_s3 = Path("c:/THz-Unified-Optimizer/docs/images/grid_integral_series3.png")
    fig_s3 = make_side_by_side_images(img_spec_s3, img_int_s3, width=245, height=120)
    if fig_s3:
        story.append(fig_s3)
        story.append(Paragraph("Рис 2. series3: спектральный анализ пропускания (слева) и интегральный метод (справа).", caption_style))
        
    # Residuals maps
    story.append(Paragraph("<b>2D-карты невязок для свободного периода P</b>", h2_style))
    img_res_amp_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_amp_356att.png")
    img_res_phase_356 = Path("c:/THz-Unified-Optimizer/docs/images/complex_residuals_phase_356att.png")
    fig_res = make_side_by_side_images(img_res_amp_356, img_res_phase_356, width=245, height=110)
    if fig_res:
        story.append(fig_res)
        story.append(Paragraph("Рис 3. Распределение амплитудных (слева) и фазовых (справа) невязок для серии 356att при свободном P.", caption_style))
        
    # Part 5
    story.append(Paragraph("5. Анализ кросс-серийной корреляции ошибок", h1_style))
    story.append(Paragraph(
        "При свободном периоде P коэффициенты глобальной взаимной корреляции остатков составили R_phi = 0.5834 (фаза) и R_A = 0.3029 (амплитуда). "
        "Тот факт, что корреляция фазовых ошибок осталась высокой (58.3%), показывает, что свободная оптимизация P и D не убирает систематические отклонения на углах скрещивания 50°--90°. "
        "Поскольку чистая апертура оправы 4 дюйма (101.6 мм) значительно превышает пучок (12 мм), дифракция на оправе исключена. Остаточные невязки вызваны ограничением динамического диапазона детектора при глубоком затухании сигнала и непараксиальностью пучка.",
        body_style
    ))
    
    img_corr_ang = Path("c:/THz-Unified-Optimizer/docs/images/correlation_angular.png")
    img_corr_scat = Path("c:/THz-Unified-Optimizer/docs/images/correlation_scatter.png")
    fig_corr = make_side_by_side_images(img_corr_ang, img_corr_scat, width=245, height=115)
    if fig_corr:
        story.append(fig_corr)
        story.append(Paragraph("Рис 4. Угловая зависимость корреляции невязок (слева) и глобальная корреляция (справа).", caption_style))
        
    # Part 6
    story.append(Paragraph("6. Итоговые выводы и калибровочные рекомендации", h1_style))
    story.append(Paragraph(
        "1. <b>Калибровка решеток:</b> При юстировке аттенюаторов нельзя рассматривать период P как жесткую паспортную константу плоской волны (15.5 мкм). Фокусировка пучка увеличивает эффективный период до 16.7--17.5 мкм. Для прецизионного моделирования необходимо использовать эффективные параметры:<br/>"
        "• 356att: P_eff = 16.729 мкм, D_eff = 4.820 мкм.<br/>"
        "• series3: P_eff = 17.500 мкм, D_eff = 4.991 мкм.<br/>"
        "2. <b>Метрологический лимит:</b> Точность Blanco-Drude модели с оптимизированным периодом достигла предела (ошибка менее 0.5 дБ для series3). Дальнейшее усложнение аналитики нецелесообразно. Для повышения точности требуется переход к 3D FDTD моделированию с учетом профиля ТГц-пучка.",
        body_style
    ))
    
    try:
        doc.build(story)
        print("Academic PDF successfully generated.")
    except PermissionError:
        pdf_path_fallback = Path("c:/THz-Unified-Optimizer/docs/artifacts/final_comprehensive_report_free_p_v2.pdf")
        doc_fallback = SimpleDocTemplate(
            str(pdf_path_fallback),
            pagesize=A4,
            rightMargin=36, leftMargin=36, topMargin=40, bottomMargin=40
        )
        doc_fallback.build(story)
        print(f"Warning: PDF file was locked. Fallback PDF generated at: {pdf_path_fallback}")

if __name__ == "__main__":
    main()
