from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from finite_slab import fastest_mode

rows = []
for eps in np.linspace(0.0, 0.99, 200):
    rows.append(fastest_mode(float(eps)))
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results" / "derived_tables" / "baseline_curves_reproduced.csv", index=False)
print(out.iloc[[0, -1]])
