import unittest
import numpy as np
from unified_optimizer.model_blanco import transmission_two_polarizers, compute_t_perp, compute_t_par

class TestModelBlanco(unittest.TestCase):
    def test_transmission_boundaries(self):
        """Проверка граничных условий модели Бланко."""
        # 1. При D -> 0 (тончайшая проволока), пропускание должно стремиться к 1.0 (ничто не мешает)
        # 2. При D -> P (проволоки сливаются), пропускание должно стремиться к 0.0 (сплошной металл)
        
        theta = 0.0
        p_over_lambda = 0.1 # квазистатический предел (длина волны много больше периода)
        
        # Граница D -> 0
        t_d_small = transmission_two_polarizers(theta, p_over_lambda, d_over_p=1e-5)
        self.assertAlmostEqual(t_d_small, 1.0, places=2)
        
        # Граница D -> P (плотное заполнение). 
        # Примечание: Аналитическая модель Бланко теряет точность при d/p -> 1
        t_d_large = transmission_two_polarizers(theta, p_over_lambda, d_over_p=0.99)
        self.assertLess(t_d_large, 0.5)
        
    def test_angle_symmetry(self):
        """Проверка симметрии T(-theta) == T(theta) и T(theta+180) == T(theta)"""
        angles = [10.0, 45.0, 75.0]
        p_over_lambda = 0.5
        d_over_p = 0.5
        
        for a_deg in angles:
            t_pos = transmission_two_polarizers(np.deg2rad(a_deg), p_over_lambda, d_over_p)
            t_neg = transmission_two_polarizers(np.deg2rad(-a_deg), p_over_lambda, d_over_p)
            t_180 = transmission_two_polarizers(np.deg2rad(a_deg + 180.0), p_over_lambda, d_over_p)
            
            self.assertAlmostEqual(t_pos, t_neg, places=5)
            self.assertAlmostEqual(t_pos, t_180, places=5)

if __name__ == '__main__':
    unittest.main()
