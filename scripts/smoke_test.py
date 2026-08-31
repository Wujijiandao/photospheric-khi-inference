from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from finite_slab import fastest_mode
from variable_density_rayleigh import solve

m = fastest_mode(0.6)
assert 0.87 < m["kappa_star"] < 0.89
assert 7.0 < m["lambda_over_d"] < 7.3
assert 0.39 < m["gamma_d_over_U0"] < 0.41

c, ga = solve(0.4446, 0.0, N=120, Lbox=12.0)
assert 0.17 < ga < 0.21
print("Smoke test passed.")
print(m)
print({"rayleigh_gamma_a_over_U0": ga, "c": c})
