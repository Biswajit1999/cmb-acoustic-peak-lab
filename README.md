# CMB Acoustic Peak Lab

Interactive cosmic microwave background power-spectrum laboratory.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of a
phenomenological model of the CMB temperature power spectrum. The app loads
`data/reference.json`, renders those published Planck-like anchor points first, then sends the
adjustable model parameters to `physicsWorker.js` so numerical work stays off the UI thread. The
goal is to build intuition for how the position, spacing and height of the acoustic peaks respond
to the underlying cosmological parameters, side by side with real benchmark values.

## Background

### The CMB as relic radiation

The cosmic microwave background is the redshifted afterglow of the last-scattering surface, the
moment (roughly 380,000 years after the Big Bang) when the universe cooled enough for electrons
and protons to combine into neutral hydrogen and photons decoupled from baryonic matter. Since
then those photons have streamed almost freely across the universe, redshifting from a plasma at
about 3000 K down to the 2.725 K blackbody we observe today. The tiny anisotropies in this
radiation, at the level of about 1 part in 100,000, are a snapshot of the density perturbations
present at decoupling, and their statistics encode the composition and geometry of the universe.

### The temperature power spectrum, C_l vs multipole l

The anisotropy field is expanded in spherical harmonics on the sky, and its statistical properties
are summarized by the angular power spectrum C_l, where the multipole moment l plays the role of
an inverse angular scale (large l corresponds to small angular scales, roughly 180 deg / l). It is
conventional to plot the rescaled quantity D_l = l(l+1)C_l / 2*pi in units of microK^2, which is
flat for scale-invariant primordial perturbations at low l and shows a series of damped
oscillations, the acoustic peaks, at higher l. This lab plots D_l against l, matching the
convention used in the Planck papers and in the reference anchors bundled with the repository.

### Baryon acoustic oscillations in the early universe

Before recombination, photons and baryons were tightly coupled by Thomson scattering, forming a
single photon-baryon fluid. Gravity from dark-matter overdensities compressed this fluid while
photon pressure resisted the compression, launching sound waves that propagated outward from each
initial perturbation at the sound speed of the fluid, c_s = c / sqrt(3(1 + R)), where
R = 3*rho_b / (4*rho_gamma) is the baryon-to-photon momentum-density ratio. At recombination the
photons decoupled and the waves froze in place, imprinted at a characteristic comoving scale, the
sound horizon r_s, equal to the distance sound could travel between the Big Bang and decoupling.
These frozen density (and corresponding temperature) modulations are the baryon acoustic
oscillations, and they show up in the CMB as a harmonic series of peaks in D_l spaced roughly
uniformly in l, at multipoles corresponding to modes that had completed a half-integer number of
compressions or rarefactions by the time of decoupling.

### How peak positions and heights constrain cosmological parameters

- **Peak spacing and the angular sound horizon.** The multipole of the first peak,
  l_1 ~ pi * D_A(z_dec) / r_s, is set by the angular size of the sound horizon on the sky, where
  D_A is the angular-diameter distance to decoupling. Because D_A depends on the universe's
  expansion history, the observed peak spacing is a sensitive probe of the spatial curvature and,
  through the distance-redshift relation, of Omega_m and H0 (equivalently h = H0 / 100 km/s/Mpc).
  A larger Omega_m or a different h shifts D_A and hence slides the whole peak pattern in l, which
  is exactly the `omegaM` and `h` sliders in this lab.
- **Peak height ratios and the baryon density.** Baryons add inertia to the photon-baryon fluid
  without adding pressure, which enhances compression peaks (odd peaks: 1st, 3rd, ...) relative to
  rarefaction peaks (even peaks: 2nd, 4th, ...). The ratio of the first to second peak heights is
  therefore one of the most direct handles on Omega_b h^2, the physical baryon density.
- **Overall tilt and the spectral index n_s.** The primordial power spectrum is close to, but not
  exactly, scale invariant. A tilt n_s != 1 changes the relative power at large versus small scales,
  tilting the whole D_l curve; this lab exposes it directly as the `ns` slider.
- **Damping tail and gravitational lensing.** Photon diffusion during recombination (Silk damping)
  exponentially suppresses power at high l, and weak gravitational lensing by intervening large-scale
  structure smooths the sharp peak features. The `lensing` slider in this lab applies an
  `exp(-lensing * l / 4800)`-type damping to caricature this smoothing so its qualitative effect on
  peak sharpness can be explored interactively.

### Connection to Planck

The Planck satellite (2009-2013 operations) produced the most precise full-sky measurement of the
CMB TT, TE and EE power spectra to date, pinning down the six-parameter LCDM model to percent-level
precision: Omega_b h^2, Omega_c h^2, H0 (or equivalently the acoustic scale theta_*), the optical
depth to reionization tau, the amplitude A_s and tilt n_s of the primordial spectrum. The reference
data bundled in `data/reference.json` is no longer a 5-point sketch -- it is the **actual 83-bin
Planck 2018 TT power spectrum**, downloaded directly from the Planck Legacy Archive
(`COM_PowerSpect_CMB-TT-binned_R3.01.txt`), with each bin's real measured `D_ell` value and
symmetrised 1-sigma uncertainty rendered as an error bar. This lets the interactive model be
visually cross-checked against the real published spectrum, all three acoustic peaks and the
damping tail included, rather than five representative landmarks.

## How It Works

This is a static, zero-build, client-side app; nothing is served except plain files.

- `index.html` — the mission-control layout: parameter sliders, an output canvas for the D_l curve,
  and a secondary heatmap canvas.
- `styles.css` — dense dark scientific-dashboard styling.
- `app.js` — owns UI state. It:
  1. Fetches `data/reference.json` on load and renders the Planck-like anchor points as labeled
     markers.
  2. Spawns `physicsWorker.js` as a Web Worker and posts the current parameter set
     (`omegaM`, `h`, `ns`, `lensing`) to it whenever a slider changes.
  3. Receives back `{ series, metrics, heatmap }` and draws the model curve and reference points
     on a shared axis (`drawSeries`), a color-mapped heatmap (`drawHeatmap`), and a metrics panel
     (peak multipole, peak power, damping ratio).
  4. Runs an independent `requestAnimationFrame` loop purely to report UI frame rate.
- `physicsWorker.js` — the numerical core, run off the main thread. For the `cmb` lab specifically,
  `cmb(p)` builds a multipole grid l in [2, 2500] and evaluates a closed-form phenomenological
  model: a smooth large-scale/Sachs-Wolfe-like rise, five Gaussian bumps standing in for the
  acoustic peaks and the damping tail (their centers scaled by a
  `(Omega_m / 0.315)^0.18 * (0.674 / h)^0.25` factor that shifts peak positions with the
  parameters), a tilt factor from `n_s`, and an exponential lensing/damping suppression at high l.
  It also produces a small 2D "acoustic response" heatmap purely for visual interest, and reports
  summary metrics (peak multipole, peak power, damping ratio, and the first-to-second peak height
  ratio, the model quantity most directly tied to the baryon density discussion below). Note that this is a fast,
  qualitatively-motivated stand-in for a Boltzmann solver (e.g. CAMB/CLASS), tuned so its shape and
  parameter sensitivities are physically sensible for teaching and exploration purposes — it does
  not solve the actual photon-baryon-dark-matter Boltzmann hierarchy. The same worker file also
  hosts unrelated toy models for other mini-labs (supernova distance moduli, FRB dispersion,
  microlensing, galaxy rotation curves, asteroseismology, weak lensing shear, spectrograph
  precision, galaxy clustering, and exoplanet transmission spectra) behind a `lab` id dispatch —
  only the `cmb` lab is wired into this app's `index.html`/`app.js`.
- `data/reference.json` — the small, auditable Planck-like reference bundle (peak locations,
  amplitudes, citation) rendered as fixed markers so the live model can be checked against them.
- `research-overlay.js` / `data/research-reference.json` / `RESEARCH_QUALITY.md` — a non-invasive
  quality/validation layer: a small mission-control panel plus benchmark anchors used by
  `scripts/validate_repository.mjs` to check repository integrity without any network access.
- `scripts/validate.js` — a no-dependency script that checks required files exist, that
  `data/reference.json` is valid JSON with finite numeric anchors, that `physicsWorker.js` parses,
  that citations are present, and that no unfinished scaffold tokens remain.

## Usage

```bash
python -m http.server 8080
```

Open `http://localhost:8080` and drag the sliders:

- `Omega_m` — total matter density parameter.
- `h` — dimensionless Hubble parameter (H0 = 100 h km/s/Mpc).
- `n_s` — primordial scalar spectral index.
- `Lensing smoothing` — strength of the high-l damping/lensing suppression applied to the model.

Watch the D_l vs l curve move relative to the fixed Planck-like reference markers, and read the
peak multipole, peak power and damping ratio off the metrics panel. Click "Reset" to return to the
Planck 2018 best-fit-like defaults (Omega_m = 0.315, h = 0.674, n_s = 0.965).

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and
absence of unfinished scaffold tokens.

## Reference Data

Planck 2018 TT power-spectrum peak locations and amplitudes, compressed for browser validation.

| Multipole l | D_l [microK^2] | Feature |
|---|---|---|
| 220 | 5740 | first acoustic peak |
| 540 | 2550 | second acoustic peak |
| 800 | 2500 | third acoustic peak |
| 1120 | 1200 | damping tail |
| 1450 | 650 | Silk damping |

## Math Appendix

Model used in `physicsWorker.js: cmb(p)`, for multipole grid l = 2..2500 (900 points):

Peak-position shift factor:

```
shift = (Omega_m / 0.315)^0.18 * (0.674 / h)^0.25
```

Spectral-tilt factor:

```
ns_factor(l) = (l / 220)^(n_s - 0.965)
```

Baseline (large-scale) term plus five Gaussian acoustic-peak bumps, each with center c_i (scaled
by `shift`), amplitude a_i and width w_i:

```
D_l = [ 850 * (l / 80)^0.2 * exp(-l / 1600)  +  sum_i a_i * exp( -0.5 * ((l - c_i) / w_i)^2 ) ]
      * ns_factor(l) * exp(-lensing * l / 4800)
```

with (c_i, a_i, w_i) = (220, 5740, 80), (540, 2550, 110), (800, 2500, 130), (1120, 1200, 170),
(1450, 650, 220) before the `shift` scaling is applied to each c_i.

For reference, the physical relations this phenomenology is meant to caricature are:

Sound speed of the photon-baryon fluid:

```
c_s = c / sqrt(3 (1 + R)),   R = 3 rho_b / (4 rho_gamma)
```

Comoving sound horizon at decoupling:

```
r_s = integral_0^{t_dec} c_s(t) dt / a(t)   =   integral_{z_dec}^{infinity} c_s(z) / H(z) dz
```

Angular position of the m-th acoustic peak (flat universe):

```
l_m ~ m * pi * D_A(z_dec) / r_s
```

D_l normalization used throughout:

```
D_l = l (l + 1) C_l / (2 pi)
```

## References

- Planck Collaboration, 2020. Planck 2018 results. VI. Cosmological parameters. Astronomy & Astrophysics, 641, A6.
- Hu, W. and Dodelson, S., 2002. Cosmic microwave background anisotropies. Annual Review of Astronomy and Astrophysics, 40, pp.171-216.

## What the Peak Positions and Ratio Actually Constrain

`peak_ell_offset_from_planck` compares the model's own first-peak location to the real Planck
2018 value (`ell ~ 220`, from the precisely measured acoustic scale `theta_* = 1.04109e-2 rad`,
Planck Collaboration, 2020, *A&A*, 641, A6, Table 2) -- at default parameters the model sits
only `~1.5` multipoles away, a direct check that the default `Omega_m`/`h` slider values are
already close to the real best-fit cosmology.

`first_second_ratio` (the first-to-second peak height ratio, `~1.79` at default parameters) is
not just a shape detail: in the real physics, **a higher baryon density `Omega_b h^2` raises
the odd (first, third, ...) peaks relative to the even (second, fourth, ...) peaks**, because
baryons add inertia to the photon-baryon plasma and load the compressions more than the
rarefactions (the same physics that makes the CMB power spectrum shape, not just its overall
amplitude, a baryon-density probe). This is precisely how Planck derives its real
`Omega_b h^2 = 0.02237 +/- 0.00015` -- from the detailed peak-height pattern, not from the peak
positions alone. Increase the `lensing` slider and watch peak heights and the damping tail
smooth out together, distinct from what changing the baryon-loading ratio would do to just the
odd/even peak heights.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
