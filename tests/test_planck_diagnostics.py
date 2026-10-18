import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from planck_diagnostics import local_peak, phenomenological_tt, residual_diagnostics, validate_points


def synthetic_points():
    ell = np.linspace(120, 330, 15)
    return [{'x': float(x), 'y': float(5800-0.16*(x-221)**2), 'y_err': 15.0} for x in ell]


def test_peak_recovers_known_vertex():
    result = local_peak(synthetic_points(), window=(120, 330), draws=400, seed=4)
    assert abs(result['ell_peak']-221) < 1e-9
    assert result['ell_16_84'][0] < 221 < result['ell_16_84'][1]


def test_peak_reproducible():
    assert local_peak(synthetic_points(), window=(120,330), draws=100,seed=12) == local_peak(synthetic_points(),window=(120,330),draws=100,seed=12)


def test_model_finite_and_positive():
    y = phenomenological_tt(np.array([2, 220, 800, 2400]))
    assert np.all(np.isfinite(y)) and np.all(y>0)


def test_model_changes_with_cosmology():
    ell=np.arange(50,2000,10)
    assert np.linalg.norm(phenomenological_tt(ell,h=.65)-phenomenological_tt(ell,h=.72))>100


def test_zero_residuals():
    points=synthetic_points()
    result=residual_diagnostics(points,model=[p['y'] for p in points],bands=((100,200),(200,400)))
    assert result['chi2']==0
    assert sum(b['bins'] for b in result['bands'])==len(points)


def test_standardized_residual():
    pts=synthetic_points()
    result=residual_diagnostics(pts,model=[p['y']+p['y_err'] for p in pts])
    assert np.isclose(result['chi2'],len(pts))


def test_duplicate_multipole_rejected():
    pts=synthetic_points();pts[1]['x']=pts[0]['x']
    with pytest.raises(ValueError):validate_points(pts)


def test_bad_error_rejected():
    pts=synthetic_points();pts[0]['y_err']=0
    with pytest.raises(ValueError):validate_points(pts)


def test_nonpeak_window_rejected():
    pts=[{'x':float(x),'y':float(x),'y_err':1.0} for x in np.arange(10,30)]
    with pytest.raises(ValueError):local_peak(pts,window=(10,30))
