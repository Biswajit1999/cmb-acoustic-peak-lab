"""Optional figure: python scripts/plot_peak_robustness.py output.png."""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.peak_robustness import audit, WINDOWS

import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
reference = json.loads((root / "data/reference.json").read_text())
points = [(p["x"], p["y"], p["y_err"]) for p in reference["points"]]
result = audit(points)
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.errorbar([p[0] for p in points], [p[1] for p in points],
            yerr=[p[2] for p in points], fmt=".", markersize=3, alpha=0.7)
for i, (fit, window) in enumerate(zip(result["peak_fits"], WINDOWS), 1):
    ax.axvspan(*window, alpha=0.06)
    ax.plot(fit["ell"], fit["height"], marker="x", markersize=8)
    ax.annotate(f"P{i}: {fit['ell']:.1f}", (fit["ell"], fit["height"]), xytext=(6, 8),
                textcoords="offset points")
ax.set(xlabel="Multipole ell", ylabel="D_ell (microK²)",
       title="Planck binned TT: local quadratic peak sensitivity")
fig.tight_layout()
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else "peak_robustness.png", dpi=170)
