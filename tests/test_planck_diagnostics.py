import math
import tempfile
import json
import unittest
from pathlib import Path
from analysis.planck_diagnostics import (toy_dl,load_planck_points,fit_amplitude,
                                         top_hat_bin_diagnostics,sensitivity_cosines,
                                         residual_by_band)

class PlanckDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.points=[(float(l),1.3*toy_dl(l),10.0) for l in range(100,2300,100)]
    def test_worker_formula_at_peak(self):
        l=220
        baseline=850*(l/80)**.2*math.exp(-l/1600)
        expected=baseline+sum(a*math.exp(-.5*((l-c)/w)**2) for c,a,w in
                              [(220,5740,80),(540,2550,110),(800,2500,130),(1120,1200,170),(1450,650,220)])
        self.assertAlmostEqual(toy_dl(l,{"omegaM":.315,"h":.674,"ns":.965,"lensing":0}),expected,8)
    def test_amplitude_exact(self):
        fit=fit_amplitude(self.points)
        self.assertAlmostEqual(fit['amplitude'],1.3,10)
        self.assertAlmostEqual(fit['chi2_diagonal'],0,8)
    def test_amplitude_with_nonuniform_errors(self):
        pts=[(x,y,5+i) for i,(x,y,s) in enumerate(self.points)]
        self.assertAlmostEqual(fit_amplitude(pts)['amplitude'],1.3,10)
    def test_amplitude_bad_input(self):
        with self.assertRaises(ValueError): fit_amplitude([(100,3,0),(200,5,1)])
        with self.assertRaises(ValueError): fit_amplitude(self.points,model_values=[0])
    def test_top_hat_integral_differs_from_center(self):
        bins=top_hat_bin_diagnostics(self.points)
        self.assertEqual(len(bins),len(self.points))
        self.assertTrue(any(abs(p['fractional_difference'])>.001 for p in bins))
    def test_top_hat_converges(self):
        a=top_hat_bin_diagnostics(self.points,steps=16)
        b=top_hat_bin_diagnostics(self.points,steps=32)
        self.assertLess(max(abs(x['top_hat_dl']-y['top_hat_dl']) for x,y in zip(a,b)),0.05)
    def test_top_hat_bad_order(self):
        with self.assertRaises(ValueError): top_hat_bin_diagnostics(list(reversed(self.points)))
    def test_omega_h_degeneracy(self):
        cos=sensitivity_cosines(self.points)['cosines']['omegaM:h']
        self.assertLess(cos,-.999999)
    def test_tilt_is_not_exactly_shift(self):
        cos=sensitivity_cosines(self.points)['cosines']['omegaM:ns']
        self.assertLess(abs(cos),.999)
    def test_residual_bands_partition(self):
        bands=residual_by_band(self.points,1.3)
        self.assertEqual(sum(x['n_bins'] for x in bands),len(self.points))
        self.assertTrue(all(abs(x['rms_normalized_residual'])<1e-9 for x in bands))
    def test_reference_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'reference.json'
            p.write_text(json.dumps({'n_points':3,'points':[{'x':a,'y':b,'y_err':c} for a,b,c in self.points[:3]]}))
            self.assertEqual(len(load_planck_points(p)),3)
            p.write_text(json.dumps({'n_points':4,'points':[{'x':a,'y':b,'y_err':c} for a,b,c in self.points[:3]]}))
            with self.assertRaises(ValueError):load_planck_points(p)
    def test_reject_invalid_model(self):
        with self.assertRaises(ValueError):toy_dl(-1)
        with self.assertRaises(ValueError):toy_dl(220,{**{'omegaM':.315,'h':.674,'ns':.965,'lensing':1},'h':-1})

if __name__=='__main__':unittest.main()
