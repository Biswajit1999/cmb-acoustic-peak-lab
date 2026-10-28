"""Run from repository root: python scripts/run_peak_robustness.py."""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.peak_robustness import audit

root = Path(__file__).resolve().parents[1]
reference = json.loads((root / "data/reference.json").read_text())
points = [(float(p["x"]), float(p["y"]), float(p["y_err"])) for p in reference["points"]]
result = audit(points)
result["input_source"] = "data/reference.json: Planck 2018 binned TT, 83 points"
print(json.dumps(result, indent=2))
