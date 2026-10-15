"""Run from repo root: python -m scripts.analyze_planck_diagnostics [reference.json]."""
import json
import sys
from pathlib import Path
from analysis.planck_diagnostics import (load_planck_points,fit_amplitude,
                                         top_hat_bin_diagnostics,sensitivity_cosines,
                                         residual_by_band)


def run(path: Path) -> dict:
    points=load_planck_points(path)
    fit=fit_amplitude(points)
    bins=top_hat_bin_diagnostics(points)
    return {'reference':str(path),'n_bins':len(points),
            'assumptions':'Toy JS model; independent symmetrized bin errors; no official Planck likelihood.',
            'diagonal_weighted_amplitude_fit':fit,
            'top_hat_approximation':{'max_abs_fractional_center_bias':max(abs(x['fractional_difference']) for x in bins),
                                     'per_bin':bins},
            'local_sensitivity':sensitivity_cosines(points),
            'residual_by_multipole_band':residual_by_band(points,fit['amplitude'])}

if __name__=='__main__':
    print(json.dumps(run(Path(sys.argv[1]) if len(sys.argv)>1 else Path('data/reference.json')),indent=2))
