import unittest
import numpy as np
from scripts.run_spatial_analysis import load_visium_sample_data, compute_morans_i

class TestSpatialAnalysis(unittest.TestCase):
    def test_load_visium_data(self):
        coords, expr, genes = load_visium_sample_data()
        self.assertEqual(len(coords), 36)
        self.assertEqual(expr.shape, (36, 2))
        
    def test_morans_i_numerical_accuracy(self):
        coords, expr, genes = load_visium_sample_data()
        score_mbp = compute_morans_i(coords, expr[:, 0])
        score_actb = compute_morans_i(coords, expr[:, 1])
        
        self.assertAlmostEqual(score_mbp, 0.4817, places=4)
        self.assertAlmostEqual(score_actb, -0.0006, places=4)

if __name__ == '__main__':
    unittest.main()
