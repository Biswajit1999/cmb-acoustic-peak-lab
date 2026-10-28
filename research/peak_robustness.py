"""Independent, local CMB TT peak diagnostics; not a cosmological likelihood.

Only diagonal, symmetrized published errors are available. In particular, the
Monte Carlo intervals below are *conditional* and omit inter-bin covariance.
"""
from __future__ import annotations

import math
import random
from statistics import median

WINDOWS = ((150.0, 290.0), (480.0, 600.0), (740.0, 890.0))


def _solve3(matrix, rhs):
    a = [list(row) + [value] for row, value in zip(matrix, rhs)]
    for i in range(3):
        pivot = max(range(i, 3), key=lambda j: abs(a[j][i]))
        if abs(a[pivot][i]) < 1e-12:
            raise ValueError("Singular quadratic fit")
        a[i], a[pivot] = a[pivot], a[i]
        v = a[i][i]
        for k in range(i, 4):
            a[i][k] /= v
        for j in range(3):
            if j != i:
                v = a[j][i]
                for k in range(i, 4):
                    a[j][k] -= v * a[i][k]
    return [a[i][3] for i in range(3)]


def peak_fit(points, window):
    """Weighted local quadratic vertex (ell, D_ell) in a declared window."""
    lo, hi = window
    selected = [p for p in points if lo <= p[0] <= hi]
    if len(selected) < 3:
        raise ValueError("At least three bins required")
    if len({p[0] for p in selected}) != len(selected):
        raise ValueError("Repeated multipoles")
    if any(not all(math.isfinite(v) for v in p) or p[2] <= 0 for p in selected):
        raise ValueError("Nonfinite data or invalid uncertainty")
    origin = (lo + hi) / 2
    # Scaled coordinates avoid a poorly conditioned raw-ell Vandermonde matrix.
    sums = [0.0] * 5
    rhs = [0.0] * 3
    for ell, dl, sigma in selected:
        t, w = (ell - origin) / 100.0, 1.0 / sigma**2
        powers = [1.0, t, t*t, t**3, t**4]
        for k in range(5):
            sums[k] += w * powers[k]
        for k in range(3):
            rhs[k] += w * dl * powers[k]
    a, b, c = _solve3([[sums[i+j] for j in range(3)] for i in range(3)], rhs)
    if c >= 0:
        raise ValueError("No concave peak in window")
    t = -b / (2*c)
    ell = origin + 100*t
    if not lo <= ell <= hi:
        raise ValueError("Fitted vertex outside window")
    return {"ell": ell, "height": a + b*t + c*t*t, "bins": len(selected)}


def _quantile(values, q):
    if not values:
        raise ValueError("Empty quantile sample")
    s = sorted(values)
    pos = (len(s)-1)*q
    i = int(pos)
    return s[i] + (s[min(i+1, len(s)-1)]-s[i])*(pos-i)


def audit(points, draws=2000, seed=20261028, windows=WINDOWS):
    """Peak localization, Gaussian error propagation, and bin deletion.

    All three windows must remain valid; failed draws are counted and omitted,
    never silently replaced. The input errors are assumed independent here.
    """
    if draws < 1:
        raise ValueError("At least one draw required")
    if len(windows) != 3:
        raise ValueError("Exactly three ordered peak windows required")
    if any(not math.isfinite(v) for p in points for v in p):
        raise ValueError("Nonfinite input")
    baseline = [peak_fit(points, w) for w in windows]
    rng = random.Random(seed)
    sampled = [[] for _ in windows]
    ratios, spacings = [], []
    rejected = 0
    for _ in range(draws):
        noisy = [(x, rng.gauss(y, s), s) for x, y, s in points]
        try:
            fits = [peak_fit(noisy, w) for w in windows]
        except ValueError:
            rejected += 1
            continue
        for dest, fit in zip(sampled, fits):
            dest.append(fit["ell"])
        ratios.append(fits[0]["height"] / fits[1]["height"])
        spacings.append(fits[2]["ell"] - fits[1]["ell"])
    deletion = []
    for w, full in zip(windows, baseline):
        selected = [i for i, p in enumerate(points) if w[0] <= p[0] <= w[1]]
        offsets, failures = [], 0
        for i in selected:
            try:
                reduced = peak_fit(points[:i] + points[i+1:], w)
                offsets.append(reduced["ell"] - full["ell"])
            except ValueError:
                failures += 1
        deletion.append({"max_abs_delta_ell": max(map(abs, offsets), default=None),
                         "valid": len(offsets), "failed": failures})
    return {
        "peak_fits": baseline,
        "gaussian_draws": draws, "accepted_draws": len(ratios), "rejected_draws": rejected,
        "peak_ell_16_50_84": [[_quantile(s, q) for q in (0.16, 0.5, 0.84)]
                               if s else None for s in sampled],
        "height_ratio_first_second_16_50_84": [_quantile(ratios, q) for q in (0.16, 0.5, 0.84)] if ratios else None,
        "third_minus_second_ell_16_50_84": [_quantile(spacings, q) for q in (0.16, 0.5, 0.84)] if spacings else None,
        "leave_one_bin_out": deletion,
        "assumptions": "Local quadratic fits; independent Gaussian symmetrized bin errors; fixed windows; no Planck covariance or cosmological inference",
    }
