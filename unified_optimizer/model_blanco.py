import numpy as np

EPS = 1e-12
C_LIGHT = 3e8
Z0 = 376.7303  # Импеданс вакуума

def estimate_deff_initial(p_phys_um: float, d_phys_um: float) -> float:
    """
    Эмпирическая оценка стартового эффективного диаметра D_eff
    на основе относительной плотности решётки (D/P):
    D_eff_init = D_phys * (1.0 - 0.85 * (D_phys / P))
    """
    d_over_p = d_phys_um / p_phys_um
    factor = np.clip(1.0 - 0.85 * d_over_p, 0.15, 0.95)
    return float(d_phys_um * factor)

def get_drude_impedance_normalized(freq_thz: float) -> complex:
    """Нормализованный поверхностный импеданс вольфрама Z_s / Z_0."""
    if freq_thz is None or freq_thz <= 0:
        return 0j
    omega = freq_thz * 2 * np.pi * 1e12
    sigma0 = 1.8e7  # См/м (Вольфрам)
    tau = 8.0e-15   # с
    sigma_omega = sigma0 / (1.0 - 1j * omega * tau)
    mu0 = 4 * np.pi * 1e-7
    Z_s = np.sqrt(1j * omega * mu0 / sigma_omega)
    return Z_s / Z0

def compute_C(m: int, p_over_lambda: float) -> complex:
    """Вычисляет параметр C_m для m-й гармоники ряда Бланко."""
    val = m**2 - p_over_lambda**2
    return np.sqrt(val + 1j * 1e-9)  # малое мнимое слагаемое против деления на 0

def safe_log(x: float) -> float:
    """Защищенный натуральный логарифм."""
    return np.log(max(x, EPS))

def compute_A1(p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """Вычисляет слагаемое A_1 для индуктивного сопротивления."""
    if d_over_p <= 0 or d_over_p >= 1:
        return 1.0
    pi_d_over_lambda = np.pi * d_over_p * p_over_lambda
    log_arg = 1.0 / (np.pi * d_over_p)
    term2 = 0.5 * (pi_d_over_lambda)**2 * (safe_log(log_arg) + 0.75)
    sum_m = 0.0
    for m in range(1, N + 1):
        sum_m += (1.0 / compute_C(m, p_over_lambda) - 1.0 / m)
    term3 = 0.5 * (pi_d_over_lambda)**2 * sum_m
    return 1.0 + term2 + term3

def compute_A2(p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """Вычисляет слагаемое A_2 для емкостного сопротивления."""
    if d_over_p <= 0 or d_over_p >= 1:
        return 1.0
    pi_d_over_lambda = np.pi * d_over_p * p_over_lambda
    log_arg = 1.0 / (np.pi * d_over_p)
    term2 = 0.5 * (pi_d_over_lambda)**2 * (11.0 / 4.0 - safe_log(log_arg))
    term3 = (1.0 / 24.0) * (pi_d_over_lambda)**2
    sum_m = 0.0
    for m in range(1, N + 1):
        sum_m += (m - 0.5 / m * (p_over_lambda)**2 - compute_C(m, p_over_lambda))
    term4 = - (pi_d_over_lambda)**2 * sum_m
    return 1.0 + term2 + term3 + term4

def compute_fa(p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """Вычисляет мнимую часть импеданса fa."""
    if d_over_p <= 0 or d_over_p >= 1:
        return 1e6
    A2 = compute_A2(p_over_lambda, d_over_p, N)
    if abs(A2) < EPS:
        return 1e6
    return 0.5 * p_over_lambda * (np.pi * d_over_p)**2 / A2

def compute_fb(p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """Вычисляет мнимую часть импеданса fb."""
    if d_over_p <= 0 or d_over_p >= 1:
        return -1e6
    A1 = compute_A1(p_over_lambda, d_over_p, N)
    A2 = compute_A2(p_over_lambda, d_over_p, N)
    if abs(A2) < EPS or abs(p_over_lambda) < EPS:
        return -1e6
    term1 = 2.0 / p_over_lambda * (1.0 / (np.pi * d_over_p))**2 * A1
    term2 = - p_over_lambda / 4.0 * (np.pi * d_over_p)**2 / A2
    return term1 + term2

def compute_fc(p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """Вычисляет комплексное сопротивление fc."""
    if d_over_p <= 0 or d_over_p >= 1:
        return 0.0
    sum_m = 0.0
    for m in range(1, N + 1):
        sum_m += (1.0 / compute_C(m, p_over_lambda) - 1.0 / m)
    log_arg = 1.0 / (np.pi * d_over_p)
    return p_over_lambda * (safe_log(log_arg) + sum_m)

def compute_fd(p_over_lambda: float, d_over_p: float) -> float:
    """Вычисляет комплексное сопротивление fd."""
    if d_over_p <= 0 or d_over_p >= 1:
        return 0.0
    return p_over_lambda * (np.pi * d_over_p)**2

def compute_t_perp(p_over_lambda: float, d_over_p: float, N: int = 15, freq_thz: float = None, use_drude: bool = True) -> complex:
    """Амплитудный коэффициент пропускания перпендикулярной поляризации."""
    fa = compute_fa(p_over_lambda, d_over_p, N)
    fb = compute_fb(p_over_lambda, d_over_p, N)
    Z_drude = get_drude_impedance_normalized(freq_thz) if use_drude else 0j
    
    Za = -1j / fa + Z_drude if abs(fa) > EPS else -1e9j
    Zb = -1j / fb + Z_drude if abs(fb) > EPS else -1e9j
    num1 = Za**2 + Za*Zb + (Za**2)*Zb
    den1 = 2.0*Za + Za**2 + Zb + Za*Zb
    Z1 = num1 / den1 if abs(den1) > EPS else 0j
    Z2 = Za / (1.0 + Za) if abs(1.0 + Za) > EPS else 1.0 + 0j
    num_t = 2.0 * Z1 * Z2
    den_t = (1.0 + Z1) * (Zb + Z2)
    return num_t / den_t if abs(den_t) > EPS else 0j

def compute_t_par(p_over_lambda: float, d_over_p: float, N: int = 15, freq_thz: float = None, use_drude: bool = True) -> complex:
    """Амплитудный коэффициент пропускания параллельной поляризации."""
    fc = compute_fc(p_over_lambda, d_over_p, N)
    fd = compute_fd(p_over_lambda, d_over_p)
    Z_drude = get_drude_impedance_normalized(freq_thz) if use_drude else 0j
    
    Zc = 1j * fc + Z_drude
    Zd = -1j * fd + Z_drude
    den34 = 1.0 + Zc + Zd
    if abs(den34) < EPS:
        return 0.0
    Z3 = (Zd + 2.0*Zc*Zd + Zd**2 + Zc) / den34
    Z4 = (Zc + Zc*Zd) / den34
    num_t = 2.0 * Z3 * Z4
    den_t = (1.0 + Z3) * (Zd + Z4) * (1.0 + Zd)
    return num_t / den_t if abs(den_t) > EPS else 0j

def transmission_two_polarizers(theta: float, p_over_lambda: float, d_over_p: float, N: int = 15) -> float:
    """
    Коэффициент пропускания по мощности для системы двух поляризаторов
    (первый развернут на угол theta в радианах). Jones matrix формализм.
    """
    t_perp = compute_t_perp(p_over_lambda, d_over_p, N)
    t_par = compute_t_par(p_over_lambda, d_over_p, N)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    R_inv = np.array([[c, s], [-s, c]])
    P = np.array([[t_perp, 0.0], [0.0, t_par]])
    M = P @ R @ P @ R_inv   # результирующая матрица системы
    E_in = np.array([1.0, 0.0])
    E_out = M @ E_in
    return float(np.clip(np.abs(E_out[0])**2 + np.abs(E_out[1])**2, 0.0, 1.0))

def transmission_complex_two_polarizers(theta: float, p_over_lambda: float, d_over_p: float, N: int = 15) -> complex:
    """
    Комплексный коэффициент пропускания для системы двух поляризаторов
    (первый развернут на угол theta в радианах). Jones matrix формализм.
    Возвращает комплексную амплитуду основной поляризации (E_out[0]).
    """
    t_perp = compute_t_perp(p_over_lambda, d_over_p, N)
    t_par = compute_t_par(p_over_lambda, d_over_p, N)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    R_inv = np.array([[c, s], [-s, c]])
    P = np.array([[t_perp, 0.0], [0.0, t_par]])
    M = P @ R @ P @ R_inv   # результирующая матрица системы
    E_in = np.array([1.0, 0.0])
    E_out = M @ E_in
    return E_out[0]


# ---------------------------------------------------------------------------
# Новая экспериментальная конфигурация: [Film P1] → [WGP(theta)] → [Film P2]
# ---------------------------------------------------------------------------

def transmission_sandwiched_wgp(theta_deg: float, p_over_lambda: float,
                                 d_over_p: float, phi1_deg: float = 0.0,
                                 phi2_deg: float = 0.0, N: int = 15,
                                 freq_thz: float = None,
                                 use_drude: bool = True) -> complex:
    """
    Комплексный коэффициент пропускания для системы:
        [Ideal Film Polarizer P1 at phi1] -> [WGP at theta] -> [Ideal Film Polarizer P2 at phi2]

    Идеальный пленочный поляризатор P_film в координатах угла phi:
        P_film(phi) = R(phi) @ [[1,0],[0,0]] @ R(-phi)

    При phi1=phi2=0 (поляризаторы вдоль оси X) результат сводится к закону Малюса:
        E_out_x = t_perp * cos^2(theta) + t_par * sin^2(theta)

    Параметры:
        theta_deg  - угол поворота WGP, градусы
        p_over_lambda - P/lambda
        d_over_p   - D/P
        phi1_deg   - ориентация входного пленочного поляризатора, градусы (обычно 0)
        phi2_deg   - ориентация выходного пленочного поляризатора, градусы (обычно 0)
        freq_thz   - частота для модели Друде
        use_drude  - использовать ли поправку Друде

    Возвращает:
        Комплексную амплитуду поля, прошедшего через систему (проекция на P2).
    """
    t_perp = compute_t_perp(p_over_lambda, d_over_p, N, freq_thz, use_drude)
    t_par  = compute_t_par (p_over_lambda, d_over_p, N, freq_thz, use_drude)

    theta = np.deg2rad(theta_deg)
    phi1  = np.deg2rad(phi1_deg)
    phi2  = np.deg2rad(phi2_deg)

    def rot(angle: float) -> np.ndarray:
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, -s], [s, c]], dtype=complex)

    # Матрица идеального плёночного поляризатора (пропускает только ось X)
    P_ideal = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)

    # Матрица WGP в собственной системе координат (ось X = перпендикулярная поляризация)
    P_wgp = np.array([[t_perp, 0.0], [0.0, t_par]], dtype=complex)

    # Входной поляризатор P1: поворачиваем в его систему, фильтруем, возвращаем
    M_P1 = rot(phi1) @ P_ideal @ rot(-phi1)

    # WGP повёрнут на theta относительно лабораторной системы
    M_WGP = rot(theta) @ P_wgp @ rot(-theta)

    # Выходной поляризатор P2
    M_P2 = rot(phi2) @ P_ideal @ rot(-phi2)

    # Полная матрица: P2 · WGP · P1
    M_total = M_P2 @ M_WGP @ M_P1

    E_in = np.array([1.0, 0.0], dtype=complex)  # Линейно поляризованный вход вдоль X
    E_out = M_total @ E_in
    # Возвращаем полную амплитуду прошедшего поля (скалярная норма для интенсивности)
    return float(np.sqrt(np.abs(E_out[0])**2 + np.abs(E_out[1])**2))


def transmission_sandwiched_wgp_complex(theta_deg: float, p_over_lambda: float,
                                         d_over_p: float, phi1_deg: float = 0.0,
                                         phi2_deg: float = 0.0, N: int = 15,
                                         freq_thz: float = None,
                                         use_drude: bool = True) -> complex:
    """
    Комплексная версия transmission_sandwiched_wgp.
    Возвращает комплексную амплитуду проекции E_out на ось P2 (вдоль phi2).
    Используется для спектрального 2D-фитинга (амплитуда + фаза).
    """
    t_perp = compute_t_perp(p_over_lambda, d_over_p, N, freq_thz, use_drude)
    t_par  = compute_t_par (p_over_lambda, d_over_p, N, freq_thz, use_drude)

    theta = np.deg2rad(theta_deg)
    phi1  = np.deg2rad(phi1_deg)
    phi2  = np.deg2rad(phi2_deg)

    def rot(angle: float) -> np.ndarray:
        c, s = np.cos(angle), np.sin(angle)
        return np.array([[c, -s], [s, c]], dtype=complex)

    P_ideal = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
    P_wgp   = np.array([[t_perp, 0.0], [0.0, t_par]], dtype=complex)

    M_P1     = rot(phi1) @ P_ideal @ rot(-phi1)
    M_WGP    = rot(theta) @ P_wgp @ rot(-theta)
    M_P2     = rot(phi2) @ P_ideal @ rot(-phi2)
    M_total  = M_P2 @ M_WGP @ M_P1

    E_in  = np.array([1.0, 0.0], dtype=complex)
    E_out = M_total @ E_in
    return E_out[0]  # Комплексная амплитуда на оси phi2
