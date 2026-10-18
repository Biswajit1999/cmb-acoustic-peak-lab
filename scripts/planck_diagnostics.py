"""Independent diagnostics for the binned Planck TT anchors in data/reference.json.

The repository's browser model is phenomenological, not a Boltzmann solver.
Results from this module must not be interpreted as cosmological posteriors.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def validate_points(points):
    """Return sorted (ell, D_ell, sigma) arrays after strict input checks."""
    a = np.array([[p['x'], p['y'], p['y_err']] for p in points], dtype=float)
    if a.ndim != 2 or a.shape[1] != 3 or len(a) < 8:
        raise ValueError('At least eight TT bins with ell, D_ell and uncertainty are needed')
    if not np.all(np.isfinite(a)) or np.any(a[:, 2] <= 0):
        raise ValueError('All TT measurements must be finite with positive errors')
    a = a[np.argsort(a[:, 0])]
    if np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError('Multipoles must be unique')
    return a.T


def local_peak(points, window=(145.0, 295.0), draws=2000, seed=0):
    """Weighted local quadratic peak and conditional noise Monte Carlo.

    The 16/84% interval assumes independent Gaussian bin errors and a locally
    quadratic spectrum. It excludes correlated systematics and model bias.
    """
    ell, power, sigma = validate_points(points)
    keep = (ell >= window[0]) & (ell <= window[1])
    if keep.sum() < 4 or draws < 20:
        raise ValueError('Peak window requires >=4 bins and >=20 draws')
    x = (ell[keep] - np.mean(ell[keep])) / 100.0
    design = np.column_stack((x*x, x, np.ones_like(x)))
    w = 1.0 / sigma[keep]
    if np.linalg.matrix_rank(design * w[:, None]) != 3:
        raise ValueError('Singular local peak fit')
    fit = np.linalg.lstsq(design*w[:, None], power[keep]*w, rcond=None)[0]
    if fit[0] >= 0:
        raise ValueError('The selected window does not contain a concave peak')
    vertex = float(np.mean(ell[keep]) - 50.0*fit[1]/fit[0])
    if not window[0] <= vertex <= window[1]:
        raise ValueError('Quadratic vertex lies outside the selected window')
    rng = np.random.default_rng(seed)
    model = design @ fit
    sample = model[:, None] + sigma[keep, None]*rng.normal(size=(keep.sum(), draws))
    simulated = np.linalg.lstsq(design*w[:, None], sample*w[:, None], rcond=None)[0]
    good = simulated[0] < 0
    simulated_vertices = np.mean(ell[keep]) - 50.0*simulated[1, good]/simulated[0, good]
    simulated_vertices = simulated_vertices[(simulated_vertices >= window[0]) & (simulated_vertices <= window[1])]
    if len(simulated_vertices) < 0.8*draws:
        raise ValueError('Peak Monte Carlo is unstable under stated errors')
    return {'ell_peak': vertex, 'ell_16_84': np.percentile(simulated_vertices, [16, 84]).tolist(),
            'n_bins': int(keep.sum()), 'accepted_draws': int(len(simulated_vertices))}


def phenomenological_tt(ell, omega_m=0.315, h=0.674, ns=0.965, lensing=0.0):
    """Python mirror of physicsWorker.js:cmb; illustrative model only."""
    ell = np.asarray(ell, dtype=float)
    if np.any(ell <= 0) or min(omega_m, h) <= 0:
        raise ValueError('Positive multipoles, omega_m and h are required')
    shift = (omega_m / 0.315)**0.18 * (0.674 / h)**0.25
    d = 850*(ell/80)**0.2*np.exp(-ell/1600)
    for center, amplitude, width in [(220,5740,80),(540,2550,110),(800,2500,130),
                                     (1120,1200,170),(1450,650,220)]:
        d += amplitude*np.exp(-0.5*((ell-center*shift)/width)**2)
    return d*(ell/220)**(ns-0.965)*np.exp(-lensing*ell/4800)


def residual_diagnostics(points, model=None, bands=((2,400),(400,1000),(1000,1800),(1800,2600))):
    """Weighted residuals and multipole-band contributions, without fitted dof claims."""
    ell, power, sigma = validate_points(points)
    pred = phenomenological_tt(ell) if model is None else np.asarray(model, dtype=float)
    if pred.shape != ell.shape or not np.all(np.isfinite(pred)):
        raise ValueError('Model must be a finite prediction for every bin')
    z = (power-pred)/sigma
    band_results = []
    for lo, hi in bands:
        mask = (ell >= lo) & (ell < hi)
        if mask.any():
            band_results.append({'ell_range': [lo, hi], 'bins': int(mask.sum()),
                                 'chi2_contribution': float(np.dot(z[mask],z[mask])),
                                 'mean_standardized_residual': float(z[mask].mean())})
    return {'n_bins': len(z), 'chi2': float(np.dot(z,z)),
            'rms_standardized_residual': float(np.sqrt(np.mean(z*z))),
            'bands': band_results, 'warning': 'Toy-model discrepancy, not a Planck likelihood'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('reference', type=Path)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    data=json.loads(args.reference.read_text())
    report={'source': data.get('dataset'), 'local_first_peak': local_peak(data['points']),
            'toy_model_residuals': residual_diagnostics(data['points'])}
    text=json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+'\n')
    else:
        print(text)

if __name__=='__main__':
    main()
