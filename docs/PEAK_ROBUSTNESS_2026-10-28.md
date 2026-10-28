# Planck TT binned-peak robustness pilot (2026-10-28)

This is an **independent descriptive audit** of `data/reference.json`, not a
replacement for a Boltzmann solver or a Planck likelihood. The browser's
phenomenological Gaussian-peak model is unchanged.

## Three numerical analyses

1. Fit weighted local quadratics to the *observed* binned `D_ell` values in
   fixed windows [150,290], [480,600], [740,890]. Vertex estimates are
   conditional on these windows and a quadratic local approximation.
2. Perturb all 83 binned values with independent normal noise using the
   **symmetrized** tabulated 1-sigma errors, refit 2,000 times with seed
   20261028, and report 16/50/84% quantiles for peak locations, first/second
   height ratio, and third-minus-second peak spacing. Invalid draws are counted.
3. Remove each bin within a window in turn and report the largest absolute
   peak-location change; this tests local-bin leverage rather than observational
   systematics or cosmic variance.

## Reproduction

```sh
python -m unittest discover -s tests -p 'test_peak_robustness.py' -v
python scripts/run_peak_robustness.py
# optional: python -m pip install matplotlib
python scripts/plot_peak_robustness.py peak_robustness.png
```

Local pre-publication read-only cross-check on the 83 committed points (1,000
rather than the script's default 2,000 perturbations): peak multipoles
**221.35, 536.04, 814.30**, peak-height ratio **2.213** (median perturbed),
third-minus-second spacing **278.4**, maximum leave-one-bin-out shifts
**2.11, 1.39, 4.67** in ell. The 1,000-draw cross-check accepted all draws.
Run the committed script for the canonical 2,000-draw output.

## Inference boundary

The actual Planck covariance matrix is **not** supplied here; bin errors
are neither independent nor a full cosmic-variance/foreground likelihood.
The perturbation quantiles are conditional computational sensitivity ranges,
**not** Planck cosmological-parameter confidence intervals. Peak spacing
cannot be directly inverted into H0, Omega_m or curvature using this local
quadratic surrogate. Do not treat the browser's parameter sliders as a
Boltzmann solver.

Source: `data/reference.json`, Planck Legacy Archive
`COM_PowerSpect_CMB-TT-binned_R3.01.txt`; Planck Collaboration (2020), A&A
641 A6. The existing reference file is not modified.
