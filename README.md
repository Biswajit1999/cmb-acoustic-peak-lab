# CMB Acoustic Peak Lab

Interactive cosmic microwave background power-spectrum laboratory.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Planck 2018 TT power-spectrum peak locations and amplitudes, compressed for browser validation.

## References

- Planck Collaboration, 2020. Planck 2018 results. VI. Cosmological parameters. Astronomy & Astrophysics, 641, A6.
- Hu, W. and Dodelson, S., 2002. Cosmic microwave background anisotropies. Annual Review of Astronomy and Astrophysics, 40, pp.171-216.
