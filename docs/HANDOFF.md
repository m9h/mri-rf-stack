# Reboot handoff: MRI RF Stack

Start here after reboot. The development work, reference results, and packaging
changes have been saved locally and published to GitHub. No running solver or
temporary download is required to use the saved results.

## What the user requested

1. Assess Hertzon plus the existing NeuroFedora COPR stack for RF simulation
   inside MRI systems.
2. Proceed with development and publish the work on GitHub.
3. Incorporate the connection to Fa-Hsuan Lin's magnetic resonance inverse imaging
   (InI), using multichannel RF head-coil sensitivities and inverse methods. The
   user confirmed Lin was the researcher they remembered from MGH and Taiwan.
4. Save the work so the machine can reboot.

The implemented release is a research prototype. Further physical validation
and realistic multichannel coil development remain; do not describe them as
completed or imply that this reproduces Lin's published acquisition results.

## Published state

| Item | Location / state at handoff |
|---|---|
| Main project | https://github.com/m9h/mri-rf-stack — branch `main` |
| Research prerelease | https://github.com/m9h/mri-rf-stack/releases/tag/v0.1.0 |
| Release source commit | `72c32fbf441fd80bbd7ac4a7ad1e248ac517488c` |
| NeuroFedora packaging | https://github.com/m9h/neurofedora/pull/2 — open, unmerged |
| Packaging branch | `feat/mri-rf-stack` |
| Packaging commit | `2895384dac54a4961ec0b9c8b141dea1cc302c44` |
| Python and Linux solver checks | https://github.com/m9h/mri-rf-stack/actions/runs/40448376813 — passed |
| Fedora 44 RPM build/install | https://github.com/m9h/mri-rf-stack/actions/runs/40448376831 — passed |

The release includes a Python wheel, Fedora 44 noarch analysis RPM and source
RPM, reference simulation ZIP, and SHA256SUMS. The handoff documentation adds a
later commit on main; the v0.1.0 source tag remains unchanged.

No COPR build was submitted. No pull request was merged. The GitHub NeuroFedora
repository is a mirror whose primary is on GitLab; the requested work was
published to GitHub.

## Local working directories

- `/Users/mhough/Workspace/mri-rf-stack`: main development repository.
- `/Users/mhough/Workspace/neurofedora-mri-rf`: isolated GitHub packaging checkout
  on `feat/mri-rf-stack`; shallow clone.
- `/Users/mhough/Workspace/neurofedora`: original GitLab checkout on
  `fix/install-blockers`. Its **pre-existing** modified `python-llvmlite.spec` and
  untracked `python-numba.spec` were left untouched. Do not discard or fold these
  unrelated edits into the MRI RF work.
- `/Users/mhough/Workspace/reboot-checkpoints/mri-rf-stack`: local recovery copy
  of the handoff, Git bundle, source snapshots, raw outputs, package artifacts,
  and the two pre-existing spec edits.

The main repo's `outputs/` and `dist/` are gitignored but persist on disk and are
included in the recovery checkpoint. The 2.3 GB `.venv` also persists on disk;
it is reproducible and is not duplicated in the recovery archive.

## Implemented and verified

- Native Hertzon eight-rung, single-feed birdcage with tuning capacitors, open
  cylindrical PEC shield, finite absorbing air boundary, and homogeneous phantom.
- Hertzon field importer preserving tetrahedron-specific values at material
  interfaces; coherent multiport accepted-power normalization.
- Algebraic B+/B− components, element-mean local SAR with polynomial integration,
  duty-cycle scaling, NPZ/JSON output and ParaView VTK export.
- JEMRIS relative complex coil-map HDF5 and XML export, explicit amplitude
  reference, phase convention, and x-fastest layout tests.
- Complex noise-whitened minimum-norm inverse, selected point-spread function,
  singular spectrum, effective degrees of freedom and noise standard deviation.
- Explicitly synthetic 32-channel encoding example. This is not a Maxwell-derived
  simulation of a 32-channel head array.
- Analysis RPM spec, CI, source/wheel builds, documented reference results.

Nineteen tests passed with Hertzon installed. Core checks passed on Python 3.11,
3.12 and 3.14; the Linux full-wave smoke solve and Fedora 44 RPM build/install
also passed. The release tag identifies the exact tested code.

## Reference results and numerical conventions

The 3 T proton benchmark frequency is 127.732435554 MHz. The final local run
stopped on Hertzon's S-parameter convergence criterion at seven passes and
16,189 tetrahedra (complex ΔS = 0.01391, threshold 0.02).

The nominal 30 pF coil is unmatched: S11 = 0.7090116305 + 0.6947668922 i;
accepted power is only 0.01460 W per 1 W incident. After scaling to 1 W accepted,
the selected phantom absorbs about 0.7524 W and has mean B+ magnitude 0.4822 µT.
These are research outputs, not validated physical performance claims.

- SI units; peak complex phasors, `exp(+iωt)`; B0 along +z.
- `Bplus = μ0(Hx+iHy)/2`, `Bminus = μ0(Hx−iHy)/2`; transmit/receive handedness
  and conjugation must match the acquisition convention.
- Accepted power: `aᴴa − (Sa)ᴴ(Sa)`, including coupling terms.
- SAR: `σ|Epeak|²/(2ρ)`. Output is element-mean local SAR, not 1 g/10 g averaging.
- Straight tet4/tet10 exports supported; curved tets and stale field caches rejected.

See `docs/validation.md` and `docs/results/report.json` for full evidence.
S-parameter convergence does not establish field/SAR convergence.

## Resume commands

The existing environment uses CPython 3.12.13. These commands work directly from
the source tree, without relying on an editable package install:

```sh
cd /Users/mhough/Workspace/mri-rf-stack
PYTHONPATH=src .venv/bin/python -m mri_rf.cli --help
MPLCONFIGDIR=outputs/.mpl PYTHONPATH=src .venv/bin/python -m pytest -q
```

Inspect existing data without recomputing:

```sh
cat outputs/convergence/analysis/report.json
```

Main output folders:

- `outputs/convergence/`: final seven-pass adaptive result, analysis, JEMRIS maps.
- `outputs/refined/`: six-pass result for comparison.
- `outputs/smoke2/`: successful one-pass local smoke solve.
- `outputs/smoke/`: earlier **failed** thin-shield-volume mesh attempt, retained
  as a diagnostic log. The corrected model uses a PEC shield surface.
- `outputs/inverse-demo/`: synthetic inverse example inputs, results, figure.
- `outputs/linux-smoke/`: downloaded successful GitHub Linux job artifacts.
- `outputs/release-artifacts/`: verified Fedora RPM and SRPM.
- `outputs/mri-rf-reference-0.1.0.zip`, `outputs/SHA256SUMS`: release copies.

To rerun, choose a new output directory; the solver intentionally rejects reuse
of a nonempty result directory:

```sh
MPLCONFIGDIR=outputs/.mpl PYTHONPATH=src .venv/bin/python -m mri_rf.cli model outputs/next-run
PYTHONPATH=src .venv/bin/python -m mri_rf.cli solve outputs/next-run/birdcage.hzn outputs/next-run/solver
PYTHONPATH=src .venv/bin/python -m mri_rf.cli analyze outputs/next-run/solver outputs/next-run/materials.json outputs/next-run/analysis
```

If rebuilding the environment, use Python 3.12, install `.[test,plot]` and
`requirements-solver.txt`, and consult the README for Linux runtime libraries.
Hertzon 0.1.25 is pinned to source
`fc82a8f64afa8ec45d44e4f8851ff39f6ffaaa70` and EMerge 3.0.0a25. Upstream requires
Python **below 3.14**; do not remove that constraint to force Fedora packaging.
`emsutil` also constrains NumPy below 2.3. The inverse implementation uses
SciPy BLAS to avoid observed spurious Apple Accelerate/NumPy complex-matmul
floating-point warnings; this path is tested.

Temporary upstream source clones and inspection files under `/private/tmp` are
re-downloadable and are not runtime dependencies of the installed environment.
The pinned Git revisions, required dependencies and exported reference inputs
are saved in the project.

## Next development work

1. Construct an independent openEMS comparison with the **same geometry, drive,
   loading, normalization and boundaries**. The stock openEMS birdcage example
   differs and is not a direct numerical reference.
2. Study field/SAR convergence, polynomial order, air-boundary distance and mesh
   sensitivity; tune and match the coil separately from checking solver plumbing.
3. Add anatomical label/mesh adapters and a realistic multichannel receive array.
4. Assemble an acquisition forward operator from complex coil sensitivities,
   reference magnetization, voxel weights and residual gradient/Fourier encoding;
   evaluate InI resolution and noise with measured/appropriate covariance.
5. Execute a full JEMRIS acquisition and reconstruction; current verification is
   export format, layout and numerical checks, not end-to-end sequence validation.
6. Audit native Hertzon/EMerge and openEMS RPM dependencies. The existing RPM is
   the analysis package only; COPR submission remains separate.

## Research connection

- Fa-Hsuan Lin et al., 2006, *Dynamic magnetic resonance inverse imaging of human
  brain function*: https://doi.org/10.1002/mrm.20997
- Lin et al., 2008, whole-brain InI using a 32-channel head coil:
  https://doi.org/10.1016/j.neuroimage.2008.04.179
- NTU profile confirming the Taiwan connection:
  https://scholars.lib.ntu.edu.tw/entities/person/fde0d05d-555b-42d7-b2b0-7f92387b0164

InI borrows inverse methods from EEG/MEG but reconstructs MR magnetization/BOLD
signals. Do not imply that RF head-coil reception measures neural currents
directly, or that the synthetic example reproduces published InI performance.
