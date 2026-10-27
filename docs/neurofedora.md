# NeuroFedora integration

The standalone analysis library is suitable for an ordinary system-Python RPM.
The initial spec is `packaging/python-mri-rf-stack.spec`; its tests run without
Hertzon. It supplies the `mri-rf` CLI, including analysis, JEMRIS export and the
inverse operator. Commands that construct or solve a Hertzon model additionally
need Hertzon in the same interpreter environment.

Existing NeuroFedora components have complementary roles:

| Packages | Integration role |
|---|---|
| 3D Slicer, SAMSEG/CHARM, FreeSurfer, ANTs, ITK | Anatomical labels, registration, geometry preparation |
| SimNIBS | Head-model preparation; its quasi-static solver is not the MRI RF engine |
| PyVista, VTK, ParaView | Inspect the exported tetrahedral field maps |
| JEMRIS | Receive exported complex relative coil maps |
| BART | Downstream reconstruction from appropriately formatted acquisition data |

## Solver packaging gate

Pinned upstream: Hertzon 0.1.25 at
`fc82a8f64afa8ec45d44e4f8851ff39f6ffaaa70`, EMerge 3.0.0a25.
Hertzon requires Python `>=3.10,<3.14`; EMerge 3 includes compiled solver kernels.
Consequently, simply dropping the upper bound and building against Fedora's
newer system Python would not be a verified packaging solution.

Its dependency chain includes EMerge, emsutil, Gmsh, NumPy, Numba/llvmlite,
PySide6, PyVista/PyVistaQt, Matplotlib, pymmg/MMG, fast-simplification, and hfsynpy.
The native kernel source build, licenses/unbundling, interpreter compatibility,
and dependency RPM availability need a Fedora build audit. The optional solver
environment is therefore pinned separately from the system analysis RPM.

Upstream also changed since the initially indexed feature matrix: 0.1.25 has
Debye/Djordjevic–Sarkar material support. This benchmark deliberately uses fixed,
explicit homogeneous material values and does not depend on those models.

openEMS has a separate native chain (fparser, TinyXML, CSXCAD, HDF5, Boost and
Python bindings). Packaging and running that solver against the same benchmark
remains the next independent-validation task. No openEMS or Hertzon binary RPM,
COPR build, or new COPR repository is claimed by this release.

## Build the core RPM

```sh
rpmbuild -ba packaging/python-mri-rf-stack.spec
```

Supply the release source archive in RPM's SOURCES directory first, and install
the spec's BuildRequires. The NeuroFedora PR carries the same spec. Its Source0
uses the public v0.1.0 source tag; it does not rely on an uncommitted local file.
