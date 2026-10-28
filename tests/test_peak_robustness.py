import math
import unittest
from research.peak_robustness import audit, peak_fit


class PeakRobustnessTests(unittest.TestCase):
    def setUp(self):
        self.windows = ((150, 290), (480, 600), (740, 890))
        self.points = []
        for center, height, lo, hi in [(220, 5800, 150, 290), (540, 2600, 480, 600), (820, 2500, 740, 890)]:
            for x in range(lo+10, hi, 20):
                self.points.append((float(x), height-0.12*(x-center)**2, 10.0))

    def test_exact_quadratic_vertex(self):
        for center, window in zip((220, 540, 820), self.windows):
            self.assertAlmostEqual(peak_fit(self.points, window)["ell"], center, places=8)

    def test_height_and_uncertainty_scale(self):
        result = audit(self.points, draws=200, seed=9, windows=self.windows)
        self.assertEqual(result["accepted_draws"], 200)
        self.assertAlmostEqual(result["peak_fits"][0]["height"], 5800, places=7)
        self.assertTrue(all(v[0] < v[2] for v in result["peak_ell_16_50_84"]))

    def test_seed_reproducibility(self):
        self.assertEqual(audit(self.points, draws=25, seed=8, windows=self.windows),
                         audit(self.points, draws=25, seed=8, windows=self.windows))

    def test_leave_one_bin_out(self):
        result = audit(self.points, draws=30, windows=self.windows)
        self.assertTrue(all(z["failed"] == 0 for z in result["leave_one_bin_out"]))
        self.assertTrue(all(z["max_abs_delta_ell"] < 1e-7 for z in result["leave_one_bin_out"]))

    def test_reject_nonconcave_and_bad_errors(self):
        with self.assertRaises(ValueError):
            peak_fit([(x, x*x, 1.0) for x in (1., 2., 3.)], (1, 3))
        with self.assertRaises(ValueError):
            peak_fit([(x, -x*x, 0.0) for x in (1., 2., 3.)], (1, 3))

    def test_reject_duplicate_ell(self):
        with self.assertRaises(ValueError):
            peak_fit([(1., 2., 1.), (1., 3., 1.), (2., 2., 1.)], (1, 3))

    def test_empty_draws(self):
        with self.assertRaises(ValueError):
            audit(self.points, draws=0)


if __name__ == "__main__":
    unittest.main()
