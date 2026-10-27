# Validation record — 2026-10-25

## Executed locally

macOS arm64, CPython 3.12.13, Hertzon 0.1.25, EMerge 3.0.0a25, Gmsh 4.15.2,
NumPy 2.2.6, SciPy 1.18.1. Hertzon source commit:
`fc82a8f64afa8ec45d44e4f8851ff39f6ffaaa70`.
See `results/environment.json` for the installed dependency versions.

Commands used for the reference run:

```sh
mri-rf model outputs/convergence --passes 12
mri-rf solve outputs/convergence/birdcage.hzn outputs/convergence/solver
mri-rf analyze outputs/convergence/solver outputs/convergence/materials.json outputs/convergence/analysis
mri-rf jemris outputs/convergence/analysis/fields.npz outputs/convergence/analysis/coil.h5 --reference-b1-ut 1
```

The same nominal model was also run with a one-pass and six-pass limit.

| Metric | 1 pass | 6 passes | 7 passes (stopped on ΔS) |
|---|---:|---:|---:|
| Tetrahedra | 3,451 | 13,224 | 16,189 |
| Selected tissue mass, kg | 4.5177 | 4.5186 | 4.5205 |
| Tissue absorbed W / 1 W accepted | 0.5495 | 0.7576 | 0.7524 |
| Mean B+ magnitude, µT / √W accepted | 0.4357 | 0.4860 | 0.4822 |
| Mean SAR, W/kg at 1 W accepted | 0.1216 | 0.1677 | 0.1664 |
| Max element-mean SAR, W/kg | 1.5370 | 2.9822 | 2.7522 |

The final complex ΔS was 0.01391, below the configured 0.02 threshold. This is
**S-parameter convergence only**. The changes in peak element-mean SAR show why
field-specific convergence must be assessed separately. This run uses low-order
elements, a coarse initial tissue mesh, and a finite absorbing air boundary.

Final S11 = 0.7090116305 + 0.6947668922 i. Only 0.01460 W is accepted per 1 W
incident, so the field normalization factor is large. The coil is deliberately
unmatched and is not claimed to be tuned to the MRI frequency. The phantom's
analytic mass is 4.5804 kg; straight-sided tetrahedra approximate its curved
surface and explain some mesh-dependent mass error.

`results/report.json` contains the complete reference metrics and all adaptive
passes. `results/birdcage.hzn`, configuration and materials reproduce the input.
The pictures show real FEM field results, not analytic/synthetic substitutions.
The separate inverse-demo picture is explicitly synthetic.

## Numerical checks

- Circularly polarized input fields select the expected algebraic B+/B− component.
- Peak-phasor SAR has the required factor 1/2.
- Coherent accepted power includes off-diagonal S-matrix coupling.
- Tetrahedral field integration matches closed-form degree-two and degree-four
  polynomial integrals and uniform-field limits.
- Accepted-power and duty-cycle scaling give the correct B1 and SAR scaling.
- Complex noise-whitened inverse agrees with the closed-form regularized solution
  and is invariant to unitary changes of channel basis.
- Singular/non-Hermitian covariance, invalid passivity, incomplete/stale output,
  incompatible geometry and material mappings fail explicitly.
- JEMRIS HDF5 sample order matches its C++ x-fastest indexing; phases are radians.

## Not yet established

No independent openEMS solve of this geometry, phantom measurement, 1 g/10 g SAR
validation, whole-body exposure, temperature prediction, or full JEMRIS sequence
execution was performed for the local reference. The InI example is a verified
linear algebra demonstration, not a reproduction of Lin et al.'s acquisition or
reported spatial/temporal resolution. No patient-safety conclusion follows from
these results.

GitHub Actions supplies separate core-Python, Linux solver-smoke, and Fedora RPM
checks. Their live statuses and logs are the authoritative record for those
platforms; the local reference above is not a claim that every hosted check ran.
