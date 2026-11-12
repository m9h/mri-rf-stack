# Hertzon–openEMS comparison

This is an independent FDTD implementation of the nominal 3 T birdcage in
`docs/results/config.json`. The first completed 6 mm run differs materially
from the original FEM reference; it is a cross-check, not a validation claim.
Mesh refinement and a higher-order FEM sensitivity check are in progress.

## First completed result

All field and absorption values use **1 W net accepted power**, peak phasors,
`exp(+iωt)`, and the same algebraic B+ definition.

| Quantity | Original Hertzon | openEMS, 6 mm |
|---|---:|---:|
| Accepted W per 1 W incident | 0.0146015 | 0.00695565 |
| Phantom mass, kg | 4.52048 | 4.63968 |
| Phantom absorption, W | 0.75237 | 0.97633 |
| Mean B+ magnitude, µT | 0.48216 | 0.55144 |
| Mean SAR, W/kg | 0.16644 | 0.21043 |

The FDTD run reached −60.10 dB energy decay after 643,795 steps. Its nominal
maximum interior cell size is 6 mm; the nonuniform grid has 627,984 Cartesian
cells (650,250 mesh nodes, which openEMS calls “FDTD cells” in its log).

On a shared 5 mm sampling grid at least 15 mm inside the phantom boundary,
the B+ magnitude NRMSE is 19.5% and complex correlation is 0.976. The complex
correlation is insensitive to a single global phase; no amplitude fit is used.
Whole-phantom means above use each solver's native integration volumes.
Peak SAR is deliberately not compared: tetrahedron means and cell-interpolated
pointwise SAR are different estimators.

## Matching and discretization

- The nominal radii, lengths, widths, gaps, eight 30 pF rung capacitors,
  homogeneous phantom properties and single 50 Ω feed are shared.
- The air box is x/y = ±320 mm and z = ±330 mm. This omits only the 0.1 µm
  OpenCASCADE tolerance in Hertzon's generated box. Hertzon's first-order
  absorbing boundary and openEMS's first-order Mur boundary are different
  numerical implementations; placing them at the same box does not establish
  identical reflection error.
- openEMS uses explicit zero-thickness Cartesian staircase sheets for the open
  shield and strips, and rasterized planar annular rings. Every rung's projected
  width is anchored at 12 mm. Mesh-dependent plane errors are recorded in the
  manifest. This is a different geometric discretization from straight FEM tets.
- Capacitors are distributed across unique axial grid edges with positive
  weights summing to 30 pF per rung. Rotated RLC boxes are not used: openEMS's
  circuit implementation snaps the untransformed box coordinates.
- The phantom is sampled at FDTD cell centres. Native material density,
  conductivity and volume are read from the raw SAR dump. Cell-interpolated
  E supplies local SAR; H is cell-interpolated in a uniformly nonmagnetic model.
- Both pulse-port and field transforms use the same one-sided DFT convention.
  Their spectral units cancel in `field/sqrt(P_acc)`. The implementation checks
  agreement between VI accepted power and the 50 Ω incident/reflected waves.
- The native PEC edge graph has four connected components: the shield, the two
  ring/rung networks, and the strip between the feed and capacitor. This checks
  that staircasing did not disconnect a rung or short the intended gaps.

## Reproduce

Install openEMS and CSXCAD C++ libraries and Python bindings separately from the
analysis package. This comparison pins the full upstream meta-repository to
`622ac90e3f37700445df424779465e6a25d3c155`, whose openEMS submodule is
`1d710c210920908a143b573f2babca05455f5479` and CSXCAD submodule is
`fdff4b551580952009c7243a29693e139d1da108`. These are development commits after
openEMS 0.37.0 and CSXCAD 0.7.0. The Python package fallback versions alone do
not identify the build; the native solver banner is saved with every run.

See the [upstream build instructions](https://docs.openems.de/en/latest/install/index.html)
and the pinned build sequence in [the Linux workflow](../.github/workflows/openems.yml).
Do not install an unrelated package named `openems` from PyPI. Use full Git
clones so native version detection retains tags. The headless build needs
CMake, a C++17 compiler, HDF5, TinyXML 1, Boost, zlib and the bundled fparser.
On macOS the upstream script builds TinyXML from checksum-verified sources.

From the project root, with the solver bindings in the same environment:

```sh
PYTHONPATH=src python -m mri_rf.openems outputs/openems-coarse \
  --config docs/results/config.json --cell-mm 6 --max-steps 2000000
PYTHONPATH=src python -m mri_rf.openems outputs/openems-fine \
  --config docs/results/config.json --cell-mm 4 --max-steps 2000000
PYTHONPATH=src python -m mri_rf.compare \
  --hertzon outputs/convergence/solver --config docs/results/config.json \
  --openems outputs/openems-coarse outputs/openems-fine --output outputs/comparison
PYTHONPATH=src python examples/plot_comparison.py outputs/comparison
```

Use `--engine gpu` for a supported native GPU backend. The default is the CPU
engine. Each solve requires a new/empty output directory; `--postprocess` reads
an existing run's saved configuration and rechecks termination. A run that hits
the step limit or is aborted is not accepted as complete. The initial wideband
pulse can excite slowly decaying modes, so the step budget is deliberately large.

`examples/check_openems_topology.py RUN_DIRECTORY` audits the exported PEC edge
graph using PyVista, available in the Hertzon environment. The comparison
evaluates FEM tet4/tet10 interpolants directly, with no averaging across element
faces, and trilinearly interpolates complex FDTD fields. Input hashes accompany
the comparison report. The Linux workflow performs a separate −40 dB coarse
smoke run; it is not the −60 dB local reference calculation.

## Interpretation limits

The coil is still unmatched, so small changes in reflection affect accepted-power
normalization strongly. Energy decay is a time-truncation check, not a mesh,
boundary or field-convergence certificate. The present differences may include
FEM order/mesh error, FDTD staircasing/material sampling, capacitor discretization
and boundary reflection. The runs do not yet separate those contributions.

No mass-averaged SAR, thermal model or measurement comparison is supplied.

Primary implementation references:

- [openEMS ports and power processing](https://docs.openems.de/en/latest/python/openEMS/ports.html)
- [Field dump interpolation conventions](https://docs.openems.de/en/latest/python/CSXCAD/CSProperties/CSPropDumpBox.html)
- [Pinned circuit discretization](https://github.com/thliebig/openEMS/blob/1d710c210920908a143b573f2babca05455f5479/FDTD/operator.cpp)
