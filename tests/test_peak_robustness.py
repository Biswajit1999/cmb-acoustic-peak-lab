"""Synthetic analytical tests for CMB local-peak inference."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import numpy as np
import pytest
from peak_robustness import quadratic_peak,peak_window_sensitivity

def test_exact_peak():
    x=np.arange(180.,271.,10.)
    y=6000-.3*(x-222)**2
    loc,err=quadratic_peak(x,y,np.full(len(x),5.))
    assert loc==pytest.approx(222,abs=1e-8)
    assert err>0

def test_error_scaling():
    x=np.arange(180.,271.,10.)
    y=6000-.3*(x-222)**2
    a=quadratic_peak(x,y,np.ones(len(x)))[1]
    b=quadratic_peak(x,y,3*np.ones(len(x)))[1]
    assert b==pytest.approx(3*a)

def test_minimum_rejected():
    x=np.arange(5.)
    with pytest.raises(ValueError):
        quadratic_peak(x,(x-2)**2,np.ones(5))

def test_window_invariance_for_exact_quadratic():
    x=np.arange(150.,301.,5.)
    y=5000-.2*(x-220)**2
    result=peak_window_sensitivity(x,y,np.full(len(x),10.),220,[25,40,65])
    assert all(abs(row[1]-220)<1e-7 for row in result)
