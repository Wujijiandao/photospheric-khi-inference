# Photospheric KHI Inference

Reproducibility package for **image-frame multimode inference from Kelvin-Helmholtz modes**, with the DKIST/MURaM photospheric detections used as an astrophysical stress test.

Accompanying revised manuscript:

**Inferring Unresolved Shear Layers from Kelvin-Helmholtz Modes: Image-Frame Multimode Identifiability and a Photospheric Stress Test**  
Yuzhan Zhang — ORCID 0009-0000-3121-7972

## What is new in v1.1.0

The central inverse no longer assumes that the layer-center velocity is independently measured.

For two **signed** modes from the same equilibrium,

- growth-rate ratios eliminate the common velocity scale;
- apparent-speed differences eliminate the unknown image-frame offset;
- two scale- and frame-free invariants recover shear thickness and density contrast under the reference model;
- the velocity scale and layer-center speed then follow algebraically.

The release also adds a scale-free global branch test. Using `r=k2/k1` and `u=kappa2/kappa_c`, the reference-slab inverse was scanned at 206,640 Jacobian points over `r=1.005–10`, `u=0.01–0.99`, and `epsilon=0.01–0.99`; the determinant remained positive throughout the sampled grid. Two independent multistart inverse surveys totaling 1,701 targets found no separated second root. These are numerical global-uniqueness tests, not an analytic theorem.

Practical conditioning still deteriorates for nearly coincident modes and near the unstable-band edges. Propagation direction must also be retained: a controlled counterexample shows that unsigned apparent-speed magnitudes can create a discrete alternative inverse branch.

A third mode overdetermines a two-parameter equilibrium family and becomes an internal closure test. Controlled smooth-profile experiments show that two modes can be structurally identifiable yet model-biased, while a withheld third mode exposes the mismatch.

## Layout

```text
src/                      Core finite-slab and smooth-profile solvers
scripts/                  Reproduction drivers
data/                     Published case values transcribed from Kuridze et al. (2026)
results/derived_tables/   Numerical outputs
figures/                  Manuscript vector figures
paper/                    Revised manuscript source and checked preview
```

## Main reproduction commands

```bash
python scripts/image_frame_multimode_analysis.py
python scripts/global_uniqueness_analysis.py --atlas --roots
python scripts/profile_mismatch_multimode.py
python scripts/muram_stress_test.py
```

The Monte Carlo and global multistart suites may take several minutes depending on hardware.

## Data provenance

The small `nature_extended_table1_*.csv` files contain values transcribed from Extended Data Table 1 of Kuridze et al. (2026), *Nature* 656, 595-601, DOI 10.1038/s41586-026-10871-3. No MURaM simulation cubes are redistributed or claimed to have been independently reanalysed.

## Versioning

- Existing archival release: v1.0.0, DOI 10.5281/zenodo.22208881
- All-versions concept DOI: 10.5281/zenodo.22208880
- This package prepares **v1.1.0**. After publishing the GitHub release, use the newly minted Zenodo version-specific DOI in the revised manuscript.

## License

Code is MIT licensed. Third-party numerical values retain their cited provenance.
