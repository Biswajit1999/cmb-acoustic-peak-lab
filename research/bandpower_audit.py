"""Planck TT binned-bandpower transforms and conditional peak-ratio errors.

This is not a substitute for the official Planck likelihood.
"""
import math

def dl_to_cl(points):
    """Convert (ell, D_ell, sigma_D) to (ell, C_ell, sigma_C)."""
    out = []
    for ell, dl, sigma in points:
        if not all(map(math.isfinite, (ell, dl, sigma))) or ell <= 0 or sigma < 0:
            raise ValueError("invalid bandpower")
        f = 2*math.pi/(ell*(ell+1))
        out.append((ell, dl*f, sigma*f))
    return out

def peak_height_ratio(peak_a, peak_b, fractional_calibration=0.0, correlation=0.0):
    """Ratio with first-order statistical error; common gain cancels exactly.

    correlation describes assumed statistical cross-bin correlation. The
    fractional_calibration argument is a common multiplicative uncertainty,
    which does not contribute to the ratio's variance.
    """
    a, sa = map(float, peak_a)
    b, sb = map(float, peak_b)
    if not all(map(math.isfinite, (a, sa, b, sb, fractional_calibration, correlation))):
        raise ValueError("nonfinite input")
    if a <= 0 or b <= 0 or sa < 0 or sb < 0 or fractional_calibration < 0 or abs(correlation) > 1:
        raise ValueError("invalid amplitude or uncertainty")
    variance = (sa/b)**2 + (a*sb/b**2)**2 - 2*correlation*a*sa*sb/b**3
    return a/b, math.sqrt(max(0.0, variance))
