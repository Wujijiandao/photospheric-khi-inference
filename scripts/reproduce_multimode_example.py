from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from finite_slab import coefficients
from multimode_inverse import invert_two_modes

# Representative truth used for the manuscript synthetic-recovery demonstration.
d_true = 10.0
u0_true = 1.5
eps_true = 0.6

# Two representative non-fastest modes, expressed directly in kappa=kd.
kappas = np.array([0.3423, 1.1638])
ks = kappas / d_true
wavelengths = 2*np.pi / ks

growth = []
vph = []
for k, kap in zip(ks, kappas):
    A, B, C = coefficients(kap, eps_true)
    D = B*B - 4*A*C
    R = -B/(2*A)
    I = np.sqrt(-D)/(2*abs(A))
    growth.append(k*u0_true*I)
    vph.append(u0_true*R)

result = invert_two_modes(wavelengths, growth, vph)
print("Synthetic truth:", {"d_km": d_true, "epsilon": eps_true, "U0_km_s": u0_true})
print("Recovered:", result)
