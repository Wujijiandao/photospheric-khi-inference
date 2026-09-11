# v1.1.0 — Image-Frame Multimode Inference

This release materially extends v1.0.0 with a general image-frame inverse for multiple Kelvin-Helmholtz modes.

## New scientific capabilities

- Derives a **two-mode image-frame inverse** that removes both the common velocity scale and the unknown shear-layer-center velocity through scale- and frame-free invariants.
- Recovers `(d, epsilon, U0, vc)` from two signed modes under a specified equilibrium family, without an external measurement of `vc`.
- Adds image-frame Jacobian and mode-pair conditioning diagnostics.
- Reparameterizes the branch problem by `r=k2/k1` and `u=kappa2/kappa_c`, making the global uniqueness test scale free.
- Adds a **206,640-point Jacobian scan** over nearly the full physical reference-slab domain. The determinant is positive at every sampled point.
- Adds **1,701 multistart inverse targets** across broad and edge-stress surveys; no separated second root is found in the tested domain.
- Explicitly maps practical ill-conditioning for nearly coincident modes and near the long-wave/cutoff boundaries.
- Demonstrates that **unsigned** apparent-speed magnitudes can generate a discrete alternative parameter branch; a common signed interface coordinate is required.
- Adds a 10,000-realization-per-setting Monte Carlo recovery test.
- Adds a three-mode closure test demonstrating the distinction between structural identifiability and profile-model adequacy.
- Adds controlled tanh, error-function, and arctangent profile-mismatch experiments.

## Photospheric application

The five published synthetic MURaM cases from Kuridze et al. (2026) remain a stress test of the fastest-mode finite-slab closure. No full MURaM cubes are redistributed or independently reanalysed.

## Reproducibility

The principal new calculations are driven by:

- `scripts/image_frame_multimode_analysis.py`
- `scripts/global_uniqueness_analysis.py`
- `scripts/profile_mismatch_multimode.py`

Derived tables, the compressed full Jacobian scan, and vector figures are included.

## Manuscript-facing update

The repository draft now includes the ApJ new-submission manuscript text prepared after the AAS80696 editorial decision. No scientific results changed after the v1.4 science freeze; the v1.5 manuscript update concerns journal-facing framing, current AAS AI/LLM disclosure, and software metadata.
