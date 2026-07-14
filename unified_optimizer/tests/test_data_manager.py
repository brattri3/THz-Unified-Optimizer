import unittest
from unified_optimizer.data_manager import parse_filename

class TestDataManager(unittest.TestCase):
    def test_parse_filename(self):
        # Обычный файл (сигнал)
        ds, angle, rep, type_ = parse_filename("series1_10deg_rep1_sig.txt")
        self.assertEqual(ds, "series1")
        self.assertEqual(angle, 10.0)
        self.assertEqual(rep, 1)
        self.assertEqual(type_, "sig")
        
        # Файл с отрицательным углом (фон)
        ds, angle, rep, type_ = parse_filename("356att_-20deg_rep1_bg.txt")
        self.assertEqual(ds, "356att")
        self.assertEqual(angle, -20.0)
        self.assertEqual(rep, 1)
        self.assertEqual(type_, "bg")
        
        # Некорректный файл
        with self.assertRaises(ValueError):
            parse_filename("wrong_format.txt")

if __name__ == '__main__':
    unittest.main()
