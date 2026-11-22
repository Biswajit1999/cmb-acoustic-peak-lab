import math
import unittest
from bandpower_audit import dl_to_cl, peak_height_ratio

class TestBandpowerAudit(unittest.TestCase):
    def test_conversion(self):
        x=dl_to_cl([(2.,6/math.pi,3/math.pi)])[0]
        self.assertAlmostEqual(x[1],2.)
        self.assertAlmostEqual(x[2],1.)
    def test_independent(self):
        r,s=peak_height_ratio((4.,.4),(2.,.2))
        self.assertAlmostEqual(r,2.)
        self.assertAlmostEqual(s,math.sqrt(.08))
    def test_common_gain(self):
        self.assertEqual(peak_height_ratio((4.,.4),(2.,.2),.3),peak_height_ratio((4.,.4),(2.,.2)))
    def test_invalid(self):
        with self.assertRaises(ValueError):
            peak_height_ratio((1.,1.),(0.,1.))
if __name__=="__main__": unittest.main()
