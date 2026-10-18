# Planck TT anchor diagnostics (independent Python analysis)

The 83 binned Planck 2018 TT measurements already pinned in `data/reference.json` are evaluated with three deliberately bounded checks:

1. **Local first-peak centroid.** Weighted quadratic in multipoles 145–295, with a seeded Gaussian-noise Monte Carlo (2000 draws). The interval is conditional on independent errors and a locally quadratic spectrum; it is **not** a cosmological-parameter uncertainty or a replacement for Planck's covariance-aware likelihood.
2. **Phenomenological-model residuals.** The Python implementation mirrors `physicsWorker.js:cmb()` at its default sliders, including its ad hoc Gaussian peaks. The resulting chi-square is an illustrative *model-discrepancy score*, not a Planck likelihood or goodness-of-fit p-value.
3. **Multipole-band residual decomposition.** Separately reports standardized-residual mean and chi-square contributions for ell 2–400, 400–1000, 1000–1800, and 1800–2600. The scores do not imply independent cosmological constraints.

```bash
python -m pip install -r requirements-research.txt
python -m pytest -q tests/test_planck_diagnostics.py
python scripts/planck_diagnostics.py data/reference.json --output results/planck_diagnostics.json
```

In a local independent run on the 83 repository bins, the first-peak vertex was ell=221.35 with conditional 16–84% noise interval [219.75, 222.95] (seed 0). The toy-model discrepancy was chi-square about 725962 for 83 bins: the model is *not* a quantitative Planck fit. Do not interpret the conditional interval as uncertainty on the true acoustic angular scale.

Input provenance: Planck 2018 binned TT file `COM_PowerSpect_CMB-TT-binned_R3.01.txt`, as documented in `data/reference.json`. This code uses the already-committed reference rather than downloading or fabricating data.
