# Photospheric KHI Inference

Reproducible theory and numerical tools for subresolution inference from photospheric Kelvin–Helmholtz (KH) modes, developed in response to the DKIST observations reported by Kuridze et al. (2026).

The repository accompanies the manuscript:

**Subresolution Inference from Photospheric Kelvin–Helmholtz Modes: Closure Limits and Multimode Identifiability**  
Yuzhan Zhang

## Scientific scope

The project asks when resolved photospheric KH modes can constrain unresolved velocity-shear layers. It contains:

- a closed finite-width slab dispersion relation and fastest-mode inversion;
- local identifiability and reference-frame degeneracy diagnostics;
- a five-case stress test using values transcribed from Extended Data Table 1 of Kuridze et al. (2026);
- smooth-profile and variable-density Rayleigh sensitivity calculations;
- profile-family sensitivity tests;
- arbitrary-mode and two-mode inversion tools;
- noisy synthetic recovery tests for multimode inference.

The main conclusion is deliberately conditional: resolved KH modes can carry subresolution information, but a universal fastest-mode wavelength-to-thickness conversion is not calibrated at event level. Profile calibration, shear-frame velocity information, or multiple modes are needed for robust inference.

## Repository layout

```text
src/                         Core numerical solvers
scripts/                     Reproduction drivers and checks
data/                        Published case values transcribed from Kuridze et al. (2026)
results/derived_tables/      Derived numerical tables
figures/                     Manuscript figures
paper/                       Manuscript source and preprint snapshot
```

## Installation

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

On macOS/Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick checks

```bash
python scripts/smoke_test.py
python scripts/reproduce_baseline.py
python scripts/reproduce_multimode_example.py
```

The committed CSV/PDF outputs correspond to the manuscript calculations. The scripts are intended to expose the numerical definitions and make the main calculations straightforward to reproduce; exact floating-point values can vary slightly with platform and dependency versions.

## Data provenance

`data/nature_extended_table1_*.csv` contains numerical values transcribed from Extended Data Table 1 of:

Kuridze, D. et al. (2026), *Nature*, **656**, 595–601, “Ubiquitous Kelvin–Helmholtz instabilities driving plasma mixing on the Sun.” DOI: 10.1038/s41586-026-10871-3.

No MURaM simulation cubes are redistributed here, and this repository does not claim a new cube-level reanalysis. The large public DKIST/MURaM data remain available from the repositories cited by Kuridze et al. (2026).

## Reproducibility boundary

The five-case MURaM comparison in the manuscript is a published-value stress test, not a blind extraction from the full radiation-MHD cubes. Smooth-profile calculations are sensitivity models and should not be interpreted as fits to the local MURaM equilibrium profiles.

## Citation

If you use this repository, cite the archived release DOI once available. GitHub also reads the included `CITATION.cff` file.

## License

Code is released under the MIT License. The manuscript and figures remain scholarly works of the author; third-party numerical values retain the provenance stated above.
