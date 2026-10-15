# Planck reference vs. browser toy-model diagnostics

Run from the repository root (no additional Python dependencies):

```sh
python -m unittest discover -s tests -p 'test_planck_diagnostics.py' -v
python -m scripts.analyze_planck_diagnostics > planck_diagnostics.json
```

The Python model reproduces the **phenomenological** `cmb()` spectrum in `physicsWorker.js`. It reads the existing `data/reference.json` with 83 Planck 2018 binned TT data points; no reference data are changed. Three analyses are provided:

1. **Weighted amplitude and residual diagnostic.** Solve analytically for the scalar amplitude minimizing `Σ[(D_obs - A D_toy)/σ]^2`, then tabulate residual RMS by four disjoint multipole ranges. The published symmetrized error bars are treated as independent for this diagnostic; **this is not a valid Planck cosmological likelihood**, and reduced χ² is not interpreted as a statistical model-selection result.
2. **Finite-bin response sensitivity.** Numerically average the toy `D_ell` over top-hat intervals whose edges are midpoints between published bin centers. Compare with center-only evaluation to measure a *model-side* discretization effect. Real Planck bandpower window functions are **not** supplied by the reference file; midpoint top-hats are explicitly approximate and must not be presented as instrument windows.
3. **Local parameter degeneracy.** Compute weighted central finite-difference derivative vectors with respect to `Omega_m`, `h`, `n_s`, and the toy lensing slider. Pairwise cosine similarities near ±1 reveal locally degenerate model responses. `Omega_m` and `h` enter only via a single power-law peak-shift combination, so their responses are exactly degenerate within this toy model; this is a limitation of the model, **not** a physical degeneracy measurement from Planck.

For a scientific cosmological parameter fit, replace this toy model with CLASS or CAMB, use the full official Planck likelihood/covariance and nuisance parameters, and propagate foreground/calibration uncertainties. These diagnostics intentionally avoid claiming a cosmological constraint.
