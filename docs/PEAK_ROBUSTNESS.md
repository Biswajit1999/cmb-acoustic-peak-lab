# Local CMB acoustic-peak sensitivity

Two independent analyses are implemented in `scripts/peak_robustness.py`:

1. A weighted quadratic fit estimates the local maximum of binned `D_ell` and propagates the coefficient covariance to the vertex. A full positive-definite bin covariance may be supplied; diagonal errors are the default.
2. Varying the fitted multipole window exposes dependence on the smooth broadband spectrum and bin selection.

These diagnostics are not Planck parameter fits, do not account for all instrumental/systematic covariance, and should not be interpreted as independent cosmological measurements. The numerical regression tests use synthetic exact quadratics, not Planck detections.

Run `python -m pytest -q tests/test_peak_robustness.py` (requires numpy and pytest).

Repository maintained by Biswajit Jana; AI-assisted coding disclosed. Scientific review is required before using these tools for inference.
