# Planck TT bandpower diagnostics (22 November 2026)

Two independent calculations: (1) transform binned D_ell to C_ell with the exact algebraic factor 2pi/[ell(ell+1)] and propagate the published symmetric bin error; (2) propagate assumed statistical cross-bin correlation into a peak-amplitude ratio, showing cancellation of a *shared* multiplicative calibration. Input: `data/reference.json` `points[]` fields x, y, y_err. Run from `research`: `python -m unittest -v test_bandpower_audit.py`.

The binned ell values are effective bin centres, not exact unbinned multipoles. A ratio of selected bins is not a fitted acoustic peak height. Statistical bin correlation is not provided by this compact JSON; do not infer it. No Planck cosmological parameter constraints are derived. Tests are synthetic. Code and tests prepared with AI assistance for review under Biswajit Jana's scientific responsibility.
