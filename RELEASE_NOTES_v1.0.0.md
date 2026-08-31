# v1.0.0 — First Archival Release

This is the first archival release of **Photospheric KHI Inference**, the reproducibility package accompanying the manuscript *Subresolution Inference from Photospheric Kelvin–Helmholtz Modes: Closure Limits and Multimode Identifiability*.

## Included

- Closed finite-width slab KHI solver and fastest-mode inversion.
- Numerical identifiability diagnostics for shear-layer thickness, velocity shear, and density contrast.
- Five-case stress test using values transcribed from Extended Data Table 1 of Kuridze et al. (2026).
- Variable-density generalized-Rayleigh solver and smooth-profile sensitivity results.
- Profile-family sensitivity calculations.
- Two-mode inversion and noisy synthetic recovery results.
- Derived tables and vector figures used in the manuscript.
- Manuscript source and a preprint snapshot.

## Scientific scope

The release supports the conclusion that resolved photospheric KH modes can encode subresolution shear-layer information, while demonstrating that a universal fastest-mode slab conversion is not calibrated at event level. Robust inference requires profile calibration, shear-frame velocity information, or multiple modes.

## Data provenance

The small tabulated DKIST/MURaM case files in `data/` are transcribed from Kuridze et al. (2026), Nature 656, 595–601, DOI: 10.1038/s41586-026-10871-3. No MURaM simulation cubes are redistributed or reanalyzed in this release.

## Reproducibility

See `README.md`, `requirements.txt`, and the scripts under `scripts/`. The committed derived tables and figures correspond to the manuscript calculations.
